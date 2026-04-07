%% DeepLabCut Automated Behavioural Analysis - Anaesthesia Emergence Phase (FBD DETECTION)
% =========================================================================
% Script associated with: "Off-target effects of DREADD ligands revealed 
% by an anaesthesia emergence paradigm in mice" (Moreno-Gomez et al., 2026)
% 
% Author: Miryam Moreno-Gómez.
% Date: April 2026
% License: MIT
%
% DESCRIPTION:
% This script identifies the First Body Movement (FBD) during the 
% emergence phase. It calculates the latency to wakefulness by detecting 
% the first coordinated movement consensus that fulfills kinematic 
% displacement and path length criteria.
%
% CITATION:
% If you use this script in your research, please cite:
% Moreno-Gómez, M. (2026). Behavioral Analysis Toolkit for FBDs. GitHub Repository.
%
% =========================================================================

clc; clear; close all;
try
    %% --- 1. CONFIGURATION AND FILE SELECTION ---
    disp('--- Step 1: Configuration and File Selection ---');
    
    % Anonymous path selection: 'pwd' starts the UI in the current folder to avoid local path leaks
    [~, path_analysis] = uigetfile('*.mat', '1/4: Select "database.mat"'); if isequal(path_analysis, 0), return; end
    [file_info, path_info] = uigetfile('*.csv', '2/4: Select TREATMENT CSV'); if isequal(file_info, 0), return; end
    [file_videos, path_videos] = uigetfile('*.xlsx', '3/4: Select VIDEOS metadata Excel'); if isequal(file_videos, 0), return; end
    path_output = uigetdir(pwd, '4/4: Select Output Folder'); if isequal(path_output, 0), return; end

    % --- DYNAMIC PARAMETERS ---
    MERGE_GAPS_UNDER_SEC = 2.0; % Gaps shorter than this (s) are merged to handle tracking flickers
    
    % --- K-Factor Analysis ---
    % The threshold is calculated as: Mean_Speed + (SD_Factor * Std_Speed)
    SD_FACTORS = [8.0]; 
    
    % Kinematic validation thresholds (adjust based on animal size/pixel scale)
    NEW_THRESHOLDS = struct('NetDisp', 50, 'PathLength', 210);
    thresholds = struct();
    for k_val = SD_FACTORS
        k_field = ['k' num2str(k_val)];
        thresholds.(k_field) = NEW_THRESHOLDS;
    end
    
    % Experimental Context
    FPS = 25; 
    LIKELIHOOD_THRESH = 0.95; % P-cut-off for DeepLabCut points
    ANALYSIS_START_SEC = 720; % Exclusion period (e.g., habituation time)
    STABLE_WINDOW_MIN_SEC = 400; % Start of baseline for noise calculation
    STABLE_WINDOW_MAX_SEC = 600; % End of baseline
    
    % Valley Splitting Logic
    % Used to separate two consecutive bouts if velocity drops below 5% of peak
    VALLEY_PEAK_RATIO = 0.05;
    MIN_VALLEY_DURATION_SEC = 3.0; 
    
    % Body Part Mapping (Indices must match your DLC config)
    BODY_PARTS_CONFIG = struct(...
        'head_centre', struct('x_row',4, 'y_row',5, 'L_row',6),...
        'body_centre', struct('x_row',16,'y_row',17,'L_row',18),...
        'tail_base',   struct('x_row',25,'y_row',26,'L_row',27));
    
    %% --- 2. MAIN ANALYSIS LOOP (K-FACTOR ITERATION) ---
    load(fullfile(path_analysis, 'database.mat'), 'database');
    if istable(database), database = table2struct(database); end
    
    % Code normalization (handles numeric or string identifiers)
    if isfield(database, 'code') && isnumeric(database(1).code)
        animalCodes_for_iteration = arrayfun(@num2str, unique([database.code]), 'UniformOutput', false);
    else
        animalCodes_for_iteration = unique({database.code});
    end

    for k_idx = 1:length(SD_FACTORS)
        SD_FACTOR = SD_FACTORS(k_idx);
        k_str_field = ['k' num2str(SD_FACTOR)];
        current_metric_thresholds = thresholds.(k_str_field);
        fprintf('\n--- Processing Analysis: K = %.1f ---\n', SD_FACTOR);
        
        % Auto-generate subfolders for organization
        k_filename_part = strrep(num2str(SD_FACTOR), '.', '_');
        plot_output_path = fullfile(path_output, ['Plots_k' k_filename_part]);
        if ~exist(plot_output_path, 'dir'), mkdir(plot_output_path); end
        trajectory_plot_path = fullfile(plot_output_path, 'Trajectories');
        if ~exist(trajectory_plot_path, 'dir'), mkdir(trajectory_plot_path); end

        final_bouts_for_this_k = {};
        for ac_idx = 1:length(animalCodes_for_iteration)
            currentAnimalCode_str = animalCodes_for_iteration{ac_idx};
            isCodeNumeric = isnumeric(database(1).code);
            if isCodeNumeric, animal_db_indices = find([database.code] == str2double(currentAnimalCode_str));
            else, animal_db_indices = find(strcmp({database.code}, currentAnimalCode_str)); end
            animalSessionDays = unique([database(animal_db_indices).day]);

            for day_idx = 1:length(animalSessionDays)
                currentSessionDay = animalSessionDays(day_idx);
                
                % Data Pre-processing: Cleans low-likelihood points and applies linear interpolation
                proc_data = process_session_data(animal_db_indices, database, currentSessionDay, BODY_PARTS_CONFIG, FPS, LIKELIHOOD_THRESH);
                if isempty(proc_data), continue; end

                % Calculate baseline movement thresholds
                movement_thresholds_vector = nan(1, 3);
                stable_period_start_frame = STABLE_WINDOW_MIN_SEC * FPS + 1;
                stable_period_end_frame = min(STABLE_WINDOW_MAX_SEC * FPS, length(proc_data.time_sec));
                if stable_period_end_frame > stable_period_start_frame
                    for bp_i = 1:3
                        stable_bp_speeds = proc_data.speeds(stable_period_start_frame:stable_period_end_frame, bp_i); 
                        valid_speeds = stable_bp_speeds(~isnan(stable_bp_speeds)); 
                        movement_thresholds_vector(bp_i) = mean(valid_speeds) + (std(valid_speeds, 0) * SD_FACTOR); 
                    end
                end
                if any(isnan(movement_thresholds_vector)), continue; end

                % Segmentation Consensus: Hierarchical approach (3 points, then 2)
                session_candidates = {};
                for points_mode = [3, 2]
                    [is_moving, movement_exceeded] = generate_is_moving_signal(proc_data.speeds, movement_thresholds_vector, 1:3, points_mode);
                    is_moving_merged = merge_short_gaps(is_moving, FPS, MERGE_GAPS_UNDER_SEC);
                    bouts = extract_bouts(is_moving_merged, proc_data, FPS, [num2str(points_mode) '-point_Merged'], movement_exceeded);
                    bouts = split_bouts_by_valley(bouts, proc_data, FPS, VALLEY_PEAK_RATIO, MIN_VALLEY_DURATION_SEC);
                    
                    for i=1:length(bouts)
                        bout = bouts{i};
                        % Validate bout against inclusion criteria
                        if bout.StartTime_sec >= ANALYSIS_START_SEC && bout.NetDisp_Mean_px > current_metric_thresholds.NetDisp && bout.PathLength_Mean_px > current_metric_thresholds.PathLength
                            session_candidates{end+1} = bout;
                        end
                    end
                    if ~isempty(session_candidates), break; end
                end
                
                if ~isempty(session_candidates)
                    % VISUALIZATION 1: VELOCITY PROFILE DIAGNOSIS
                    fig1 = figure('Visible', 'off', 'Position', [100 100 1200 600]);
                    hold on;
                    colors = [1 0 0; 0 0.7 0; 0 0 1]; 
                    body_parts = fieldnames(BODY_PARTS_CONFIG);
                    for bp_i = 1:3, plot(proc_data.time_sec, proc_data.speeds(:, bp_i), 'Color', [colors(bp_i, :) 0.4], 'LineWidth', 1.5, 'DisplayName', ['Vel. ' body_parts{bp_i}]); end
                    for bp_i = 1:3, plot([proc_data.time_sec(1), proc_data.time_sec(end)], [movement_thresholds_vector(bp_i), movement_thresholds_vector(bp_i)], '--', 'Color', colors(bp_i, :), 'LineWidth', 2, 'DisplayName', ['Threshold ' body_parts{bp_i}]); end
                    y_limits = get(gca, 'YLim');
                    for b_plot_idx = 1:length(session_candidates), bout_to_plot = session_candidates{b_plot_idx}; startTime = bout_to_plot.StartTime_sec; endTime = startTime + bout_to_plot.Duration_sec; patch([startTime endTime endTime startTime], [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], 'k', 'FaceAlpha', 0.15, 'EdgeColor', 'none'); end
                    p_dummy = patch(NaN, NaN, 'k', 'FaceAlpha', 0.15, 'EdgeColor', 'none'); set(get(get(p_dummy,'Annotation'),'LegendInformation'),'IconDisplayStyle','on'); set(p_dummy, 'DisplayName', 'Detected Bouts');
                    title_str = sprintf('Session Profile - Animal: %s / Day: %d / k = %.1f', currentAnimalCode_str, currentSessionDay, SD_FACTOR);
                    title(title_str, 'Interpreter', 'none'); xlabel('Time (s)'); ylabel('Velocity (px/s)'); legend('show', 'Location', 'northeast'); grid on; hold off;
                    plot_filename_1 = fullfile(plot_output_path, sprintf('VelocityProfile_Animal%s_Day%d.svg', currentAnimalCode_str, currentSessionDay));
                    try, saveas(fig1, plot_filename_1); catch, end;
                    close(fig1);

                    % CONSISTENT AXIS CALCULATION (for standardized trajectory plotting)
                    session_min_x = Inf; session_max_x = -Inf;
                    session_min_y = Inf; session_max_y = -Inf;
                    for b_idx = 1:length(session_candidates)
                        bout = session_candidates{b_idx};
                        start_frame = round(bout.StartTime_sec * FPS) + 1;
                        end_frame = start_frame + round(bout.Duration_sec * FPS) - 1;
                        end_frame = min(end_frame, size(proc_data.coords, 1));
                        if start_frame >= end_frame, continue; end
                        bout_coords = proc_data.coords(start_frame:end_frame, :);
                        all_x_coords = bout_coords(:, 1:2:end);
                        all_y_coords = bout_coords(:, 2:2:end);
                        session_min_x = min(session_min_x, min(all_x_coords(:)));
                        session_max_x = max(session_max_x, max(all_x_coords(:)));
                        session_min_y = min(session_min_y, min(all_y_coords(:)));
                        session_max_y = max(session_max_y, max(all_y_coords(:)));
                    end
                    
                    x_range = session_max_x - session_min_x; if x_range == 0, x_range = 10; end
                    y_range = session_max_y - session_min_y; if y_range == 0, y_range = 10; end
                    x_padding = x_range * 0.1; y_padding = y_range * 0.1;
                    axis_limits = [session_min_x - x_padding, session_max_x + x_padding, session_min_y - y_padding, session_max_y + y_padding];
                    
                    % VISUALIZATION 2: SPATIAL TRAJECTORIES
                    for b_idx = 1:length(session_candidates)
                        generate_bout_trajectory_plot(session_candidates{b_idx}, proc_data, trajectory_plot_path, FPS, BODY_PARTS_CONFIG, currentAnimalCode_str, currentSessionDay, axis_limits);
                    end
                end

                if isempty(session_candidates), continue; end
                
                % SELECTION LOGIC: Choosing the most representative bout for final reporting
                if length(session_candidates) == 1, final_bout_for_session = session_candidates{1};
                else
                    products = cellfun(@(c) c.NetDisp_Mean_px * c.PathLength_Mean_px, session_candidates);
                    [~, max_idx] = max(products);
                    final_bout_for_session = session_candidates{max_idx};
                end
                
                final_bout_for_session.AnimalCode = str2double(currentAnimalCode_str);
                final_bout_for_session.SessionDay = currentSessionDay;
                final_bouts_for_this_k{end+1} = final_bout_for_session;
            end
        end
        
        % FINAL DATA EXPORT
        fprintf('--- Generating output table (K=%.1f) ---\n', SD_FACTOR);
        if ~isempty(final_bouts_for_this_k)
            T_k = struct2table([final_bouts_for_this_k{:}]);
            T_output = create_output_table(T_k, path_info, file_info, path_videos, file_videos);
            output_filename = fullfile(path_output, ['Summary_Results_k' k_filename_part '.xlsx']);
            writetable(T_output, output_filename, 'Sheet', 'Kinematic_Data');
            fprintf('Results successfully saved: %s\n', output_filename);
        else
            fprintf('No bouts detected for K = %.1f.\n', SD_FACTOR);
        end
    end
    fprintf('\n--- PROCESS COMPLETED ---\n');
catch ME, fprintf(2, '\nCRITICAL ERROR: %s, Line %d\n', ME.message, ME.stack(1).line); end

%% --- AUXILIARY FUNCTIONS ---

function T_output = create_output_table(T_k, path_info, file_info, path_videos, file_videos)
    % MERGES KINEMATIC RESULTS WITH TREATMENT AND VIDEO METADATA
    if isempty(T_k) || height(T_k) == 0, T_output = table(); return; end
    T_info_tratamiento = readtable(fullfile(path_info, file_info), 'Delimiter', ';');
    var_names = T_info_tratamiento.Properties.VariableNames;
    var_names{strcmp(var_names, 'Animal')} = 'AnimalCode';
    var_names{strcmp(var_names, 'Grupo')} = 'Treatment';
    T_info_tratamiento.Properties.VariableNames = var_names;
    if iscell(T_info_tratamiento.AnimalCode), T_info_tratamiento.AnimalCode = cellfun(@(x) str2double(strrep(x, 'Bc', '')), T_info_tratamiento.AnimalCode); end
    if isnumeric(T_info_tratamiento.Week_day), T_info_tratamiento.Week_day = cellstr(num2str(T_info_tratamiento.Week_day)); end
    partes = split(T_info_tratamiento.Week_day, '.');
    T_info_tratamiento.SessionDay = str2double(partes(:, 1));
    mapa_tratamiento = unique(T_info_tratamiento(:, {'AnimalCode', 'SessionDay', 'Treatment'}), 'rows');
    T_info_videos = readtable(fullfile(path_videos, file_videos));
    var_names = T_info_videos.Properties.VariableNames;
    var_names{strcmp(var_names, 'code')} = 'AnimalCode';
    var_names{strcmp(var_names, 'day')} = 'SessionDay';
    T_info_videos.Properties.VariableNames = var_names;
    tokens = regexp(T_info_videos.xy, '(\d+)DLC', 'tokens', 'once');
    T_info_videos.VideoID = cellfun(@(c) c{1}, tokens, 'UniformOutput', false);
    if iscell(T_info_videos.AnimalCode), T_info_videos.AnimalCode = cellfun(@(x) str2double(strrep(x, 'Bc', '')), T_info_videos.AnimalCode); end
    [claves_agrupadas, ~, ic] = unique(T_info_videos(:, {'AnimalCode', 'SessionDay'}), 'rows');
    lista_videos_agrupados = cell(height(claves_agrupadas), 1);
    for i = 1:height(claves_agrupadas)
        lista_videos_agrupados{i} = strjoin(T_info_videos.VideoID(ic == i), ', ');
    end
    claves_agrupadas.Videos = lista_videos_agrupados;
    T_k_enriquecida = outerjoin(T_k, mapa_tratamiento, 'Keys', {'AnimalCode', 'SessionDay'}, 'Type', 'left', 'MergeKeys', true);
    T_k_enriquecida = outerjoin(T_k_enriquecida, claves_agrupadas, 'Keys', {'AnimalCode', 'SessionDay'}, 'Type', 'left', 'MergeKeys', true);
    final_cols = { 'AnimalCode', 'SessionDay', 'Treatment', 'Videos', 'StartTime_sec', 'Duration_sec', 'NetDisp_Mean_px', 'PathLength_Mean_px', 'MaxSpeed_px_s', 'MinSpeed_px_s', 'AUC_px', 'DetectionType', 'InvolvedPoints' };
    T_output = T_k_enriquecida(:, final_cols);
    T_output.Properties.VariableNames = { 'Animal', 'Day', 'Treatment', 'Videos', 'Timestamp_s', 'Duration_s', 'NetDisplacement_px', 'PathLength_px', 'MaxSpeed_px_s', 'MinSpeed_px_s', 'AUC_Velocity', 'DetectionMethod', 'ActivePoints' };
end

function proc_data = process_session_data(db_indices, database, day, config, fps, L_thresh)
    % PERFORMS DATA CLEANING AND LINEAR INTERPOLATION OF TRACKING GAPS
    proc_data = []; 
    session_video_db_indices = db_indices([database(db_indices).day] == day);
    if ~isempty(session_video_db_indices)
        bodyPartNames = fieldnames(config);
        numBodyParts = length(bodyPartNames);
        all_X = cell(1, numBodyParts); all_Y = cell(1, numBodyParts); all_L = cell(1, numBodyParts);
        current_offset = 0;
        for i = 1:length(session_video_db_indices)
            h5_path = strtrim(database(session_video_db_indices(i)).xy);
            if ~exist(h5_path, 'file'), continue; end
            try
                dlc_data = h5read(h5_path, '/df_with_missing/table/');
                num_frames = size(dlc_data.values_block_0, 2);
                for bp_idx = 1:numBodyParts
                    bp_name = bodyPartNames{bp_idx}; bp_c = config.(bp_name);
                    all_X{bp_idx} = [all_X{bp_idx}; dlc_data.values_block_0(bp_c.x_row, :)'];
                    all_Y{bp_idx} = [all_Y{bp_idx}; dlc_data.values_block_0(bp_c.y_row, :)'];
                    all_L{bp_idx} = [all_L{bp_idx}; dlc_data.values_block_0(bp_c.L_row, :)'];
                end
                current_offset = current_offset + num_frames;
            catch
                warning('Error reading H5: %s', h5_path);
            end
        end
        if current_offset > 0
            proc_data.coords = NaN(current_offset, numBodyParts * 2);
            proc_data.speeds = NaN(current_offset, numBodyParts);
            for bp_idx = 1:numBodyParts
                x = all_X{bp_idx}; y = all_Y{bp_idx}; L = all_L{bp_idx};
                invalid_idx = L < L_thresh;
                x(invalid_idx) = NaN; y(invalid_idx) = NaN;
                nan_indices_x = find(isnan(x)); valid_indices_x = find(~isnan(x));
                if length(valid_indices_x) >= 2, interp_x = x; interp_x(nan_indices_x) = interp1(valid_indices_x, x(valid_indices_x), nan_indices_x, 'linear', 'extrap'); interp_x = fillmissing(interp_x, 'nearest');
                else, interp_x = fillmissing(x, 'nearest'); end
                nan_indices_y = find(isnan(y)); valid_indices_y = find(~isnan(y));
                if length(valid_indices_y) >= 2, interp_y = y; interp_y(nan_indices_y) = interp1(valid_indices_y, y(valid_indices_y), nan_indices_y, 'linear', 'extrap'); interp_y = fillmissing(interp_y, 'nearest');
                else, interp_y = fillmissing(y, 'nearest'); end
                proc_data.coords(:, (bp_idx*2)-1) = interp_x;
                proc_data.coords(:, bp_idx*2) = interp_y;
                proc_data.speeds(:, bp_idx) = [NaN; sqrt(diff(interp_x).^2 + diff(interp_y).^2) * fps];
            end
            proc_data.time_sec = (0:current_offset-1)' / fps;
        end
    end
end

function [is_moving, movement_exceeded_per_bp] = generate_is_moving_signal(speeds, thresholds, main_bp_indices, min_bps_moving)
    % LOGICAL SIGNAL GENERATOR FOR MOVEMENT CONSENSUS
    main_speeds = speeds(:, main_bp_indices);
    movement_exceeded_per_bp = main_speeds > thresholds;
    num_bps_moving = sum(movement_exceeded_per_bp, 2);
    is_moving = num_bps_moving >= min_bps_moving;
end

function all_bouts = extract_bouts(is_moving, proc_data, fps, detection_type, movement_exceeded_per_bp)
    % EXTRACTS DISCRETE SEGMENTS FROM CONTINUOUS BINARY SIGNALS
    all_bouts = {}; bp_names_map = {'H', 'B', 'T'};
    movement_starts = find(diff([false; is_moving]) == 1);
    movement_ends = find(diff([is_moving; false]) == -1);
    if isempty(movement_starts) || isempty(movement_ends), return; end
    if length(movement_starts) > length(movement_ends), movement_starts(end) = []; end
    if length(movement_ends) > length(movement_starts), movement_ends(1) = []; end
    if isempty(movement_starts) || isempty(movement_ends), return; end
    if movement_starts(1) > movement_ends(1), movement_ends(1) = []; end
    if isempty(movement_starts) || isempty(movement_ends), return; end
    if movement_starts(end) > movement_ends(end), movement_starts(end) = []; end
    for i = 1:length(movement_starts)
        start_idx = movement_starts(i); end_idx = movement_ends(i);
        if end_idx <= start_idx, continue; end
        s = struct(); s.StartTime_sec = proc_data.time_sec(start_idx);
        s.Duration_sec = (end_idx - start_idx + 1) / fps;
        if isempty(movement_exceeded_per_bp), involved_points_indices = 1:3; else, bout_movement_matrix = movement_exceeded_per_bp(start_idx:end_idx, :); involved_points_indices = find(any(bout_movement_matrix, 1)); end
        if isempty(involved_points_indices), involved_points_indices = 1:3; end
        bout_coords = proc_data.coords(start_idx:end_idx, :);
        net_disps = zeros(1, length(involved_points_indices));
        for j = 1:length(involved_points_indices), k = involved_points_indices(j); net_disps(j) = sqrt((bout_coords(end, (k*2)-1) - bout_coords(1, (k*2)-1))^2 + (bout_coords(end, k*2) - bout_coords(1, k*2))^2); end
        s.NetDisp_Mean_px = mean(net_disps);
        path_lengths = zeros(1, length(involved_points_indices));
        for j = 1:length(involved_points_indices), k = involved_points_indices(j); path_lengths(j) = sum(sqrt(diff(bout_coords(:, (k*2)-1)).^2 + diff(bout_coords(:, k*2)).^2), 'omitnan'); end
        s.PathLength_Mean_px = mean(path_lengths);
        bout_speeds = proc_data.speeds(start_idx:end_idx, :);
        mean_bout_speed_profile = mean(bout_speeds(:, involved_points_indices), 2, 'omitnan');
        s.MaxSpeed_px_s = max(mean_bout_speed_profile);
        s.MinSpeed_px_s = min(mean_bout_speed_profile);
        s.AUC_px = trapz(mean_bout_speed_profile) * (1 / fps);
        s.DetectionType = detection_type;
        s.InvolvedPoints = strjoin(bp_names_map(involved_points_indices), ',');
        all_bouts{end+1} = s;
    end
end

function final_bouts = split_bouts_by_valley(bouts_to_check, proc_data, fps, valley_peak_ratio, min_valley_sec)
    % REFINEMENT STEP: SPLITS COMPOUND BOUTS SEPARATED BY SHORT PAUSES
    if isempty(bouts_to_check), final_bouts = {}; return; end
    final_bouts = {}; min_valley_frames = min_valley_sec * fps;
    for i = 1:length(bouts_to_check)
        current_bout = bouts_to_check{i};
        start_frame = round(current_bout.StartTime_sec * fps) + 1;
        end_frame = start_frame + round(current_bout.Duration_sec * fps) - 1;
        end_frame = min(end_frame, size(proc_data.speeds, 1));
        bout_mean_speed = mean(proc_data.speeds(start_frame:end_frame, :), 2, 'omitnan');
        peak_speed = max(bout_mean_speed);
        if peak_speed == 0, final_bouts{end+1} = current_bout; continue; end
        valley_threshold = peak_speed * valley_peak_ratio;
        is_in_valley = bout_mean_speed < valley_threshold;
        valley_starts = find(diff([false; is_in_valley]) == 1);
        valley_ends = find(diff([is_in_valley; false]) == -1);
        split_points_frames = [];
        if ~isempty(valley_starts)
            for j = 1:length(valley_starts)
                if valley_ends(j) - valley_starts(j) + 1 >= min_valley_frames, split_points_frames = [split_points_frames; valley_starts(j), valley_ends(j)]; end
            end
        end
        if isempty(split_points_frames), final_bouts{end+1} = current_bout;
        else
            last_cut_frame = 1;
            for j = 1:size(split_points_frames, 1)
                new_bout_start_abs = start_frame + last_cut_frame - 1;
                new_bout_end_abs = start_frame + split_points_frames(j, 1) - 2;
                if new_bout_end_abs > new_bout_start_abs
                   is_moving_segment = false(size(proc_data.speeds, 1), 1); is_moving_segment(new_bout_start_abs:new_bout_end_abs) = true;
                   recalculated_bout = extract_bouts(is_moving_segment, proc_data, fps, current_bout.DetectionType, []);
                   if ~isempty(recalculated_bout), final_bouts = [final_bouts, recalculated_bout]; end
                end
                last_cut_frame = split_points_frames(j, 2) + 1;
            end
            new_bout_start_abs = start_frame + last_cut_frame - 1;
            new_bout_end_abs = end_frame;
            if new_bout_end_abs > new_bout_start_abs
               is_moving_segment = false(size(proc_data.speeds, 1), 1); is_moving_segment(new_bout_start_abs:new_bout_end_abs) = true;
               recalculated_bout = extract_bouts(is_moving_segment, proc_data, fps, current_bout.DetectionType, []);
               if ~isempty(recalculated_bout), final_bouts = [final_bouts, recalculated_bout]; end
            end
        end
    end
end

function merged_signal = merge_short_gaps(is_moving, fps, max_gap_sec)
    % TEMPORAL FILTERING: CONNECTS SEGMENTS BROKEN BY TRACKING NOISE
    max_gap_frames = round(max_gap_sec * fps);
    merged_signal = is_moving;
    pause_starts = find(diff([true; ~is_moving]) == 1);
    pause_ends = find(diff([~is_moving; true]) == -1);
    if isempty(pause_starts) || isempty(pause_ends), return; end
    if pause_ends(1) < pause_starts(1), pause_ends(1) = []; end
    if isempty(pause_ends), return; end
    if pause_starts(end) > pause_ends(end), pause_starts(end) = []; end
    for i = 1:length(pause_starts)
        if pause_ends(i) - pause_starts(i) + 1 <= max_gap_frames, merged_signal(pause_starts(i):pause_ends(i)) = true; end
    end
end

function generate_bout_trajectory_plot(bout, proc_data, output_path, fps, config, animal_code, session_day, axis_limits)
    % GENERATES HIGH-RESOLUTION SPATIAL PLOTS FOR INDIVIDUAL BOUTS
    start_frame = round(bout.StartTime_sec * fps) + 1;
    end_frame = start_frame + round(bout.Duration_sec * fps) - 1;
    end_frame = min(end_frame, size(proc_data.coords, 1));
    if start_frame >= end_frame, return; end
    bout_coords = proc_data.coords(start_frame:end_frame, :);

    involved_points_str = strsplit(bout.InvolvedPoints, ',');
    x_coords_to_avg = []; y_coords_to_avg = [];
    for i = 1:length(involved_points_str)
        point_char = involved_points_str{i}; bp_idx = 0;
        if contains(point_char, 'H'), bp_idx = 1; end
        if contains(point_char, 'B'), bp_idx = 2; end
        if contains(point_char, 'T'), bp_idx = 3; end
        if bp_idx > 0
            x_coords_to_avg = [x_coords_to_avg, bout_coords(:, (bp_idx*2)-1)];
            y_coords_to_avg = [y_coords_to_avg, bout_coords(:, bp_idx*2)];
        end
    end
    if isempty(x_coords_to_avg), mean_x = mean(bout_coords(:, 1:2:end), 2, 'omitnan'); mean_y = mean(bout_coords(:, 2:2:end), 2, 'omitnan');
    else, mean_x = mean(x_coords_to_avg, 2, 'omitnan'); mean_y = mean(y_coords_to_avg, 2, 'omitnan'); end
    
    fig = figure('Visible', 'off', 'Position', [100 100 800 800]);
    hold on; grid on; axis equal;
    colors = [1 0 0; 0 0.7 0; 0 0 1]; body_parts = fieldnames(config);
    
    for bp_i = 1:3
        x = bout_coords(:, (bp_i*2)-1); y = bout_coords(:, bp_i*2);
        plot(x, y, 'Color', [colors(bp_i, :) 0.3], 'LineWidth', 1, 'DisplayName', body_parts{bp_i});
        plot([x(1) x(end)], [y(1) y(end)], '--', 'Color', [colors(bp_i, :) 0.5], 'HandleVisibility','off');
        scatter(x(1), y(1), 30, 'o', 'filled', 'MarkerFaceColor', colors(bp_i,:), 'HandleVisibility','off');
        scatter(x(end), y(end), 30, 'x', 'LineWidth', 1.5, 'MarkerEdgeColor', colors(bp_i,:), 'HandleVisibility','off');
    end
    
    plot(mean_x, mean_y, 'k-', 'LineWidth', 2.5, 'DisplayName', 'Active Points (Mean)');
    plot([mean_x(1) mean_x(end)], [mean_y(1) mean_y(end)], 'k--','HandleVisibility','off');
    scatter(mean_x(1), mean_y(1), 100, 'o', 'filled', 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'k', 'DisplayName', 'Bout Start');
    scatter(mean_x(end), mean_y(end), 100, 'x', 'LineWidth', 3, 'MarkerEdgeColor', 'r', 'DisplayName', 'Bout End');
    
    title_str = sprintf('Trajectory Bout (%.2fs) - Animal %s, Day %d\nDisp: %.2f px, Path: %.2f px', bout.StartTime_sec, animal_code, session_day, bout.NetDisp_Mean_px, bout.PathLength_Mean_px);
    title(title_str, 'Interpreter', 'none'); xlabel('X Coordinate (px)'); ylabel('Y Coordinate (px)');
    legend('show', 'Location', 'best');
    axis(axis_limits);
    set(gca, 'YDir','reverse'); % Typical for video-based coordinate systems
    hold off;
    
    filename = fullfile(output_path, sprintf('Trajectory_Animal%s_Day%d_Bout%.2fs.svg', animal_code, session_day, bout.StartTime_sec));
    try, saveas(fig, filename); catch, end;
    close(fig);
end