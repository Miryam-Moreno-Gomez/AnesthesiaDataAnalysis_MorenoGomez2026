%% DeepLabCut Automated Behavioural Analysis - Anaesthesia Induction Phase (LBM DETECTION)
% =========================================================================
% Script associated with: "Off-target effects of DREADD ligands revealed 
% by an anaesthesia emergence paradigm in mice" (Moreno-Gomez et al., 2026)
% 
% Author: Miryam Moreno-Gómez.
% Date: April 2026
% License: MIT
% 
% OBJECTIVE:
% Identifies the 'Last Body Movement' (LBM) during the anaesthesia 
% induction phase (transition into isoflurane-induced sleep).
%
% CITATION:
% If you use this script in your research, please cite:
% Moreno-Gómez, M. (2026). Behavioral Analysis Toolkit for FBDs. GitHub Repository.
%

clc; clear; close all;

%% --- 1. CONFIGURATION AND ANALYSIS PARAMETERS ---
FPS = 25;                       % Frames per second [cite: 78]
LIKELIHOOD_THRESH = 0.95;       % DeepLabCut confidence threshold [cite: 82]
SESSION_TRIM_SEC = 15;          % Seconds to trim from start/end of session
ISO_ADMIN_START_SEC = 120;      % Start of Isoflurane administration (2 min baseline) [cite: 68]
LBM_WINDOW_END_SEC = 540;       % End of the induction search window (9 min)

% Body Parts Configuration (Standard 10-label model)
BODY_PARTS_CONFIG = struct(...
    'snout',           struct('name', 'snout',           'x_row', 1,  'y_row', 2,  'L_row', 3), ...
    'head_centre',     struct('name', 'head_centre',     'x_row', 4,  'y_row', 5,  'L_row', 6), ...
    'left_ear',        struct('name', 'left_ear',        'x_row', 7,  'y_row', 8,  'L_row', 9), ...
    'right_ear',       struct('name', 'right_ear',       'x_row', 10, 'y_row', 11, 'L_row', 12), ...
    'neck',            struct('name', 'neck',            'x_row', 13, 'y_row', 14, 'L_row', 15), ...
    'body_centre',     struct('name', 'body_centre',     'x_row', 16, 'y_row', 17, 'L_row', 18), ...
    'left_body_side',  struct('name', 'left_body_side',  'x_row', 19, 'y_row', 20, 'L_row', 21), ... 
    'right_body_side', struct('name', 'right_body_side', 'x_row', 22, 'y_row', 23, 'L_row', 24), ... 
    'tail_base',       struct('name', 'tail_base',       'x_row', 25, 'y_row', 26, 'L_row', 27), ... 
    'tail_tip',        struct('name', 'tail_tip',        'x_row', 28, 'y_row', 29, 'L_row', 30)  ... 
);
bodyPartNames = fieldnames(BODY_PARTS_CONFIG);

% Detection Logic Indices
idx_hc = find(strcmp('head_centre', bodyPartNames));
idx_bc = find(strcmp('body_centre', bodyPartNames));
idx_tb = find(strcmp('tail_base', bodyPartNames));

% Threshold Config (T_Mean_3_0SD)
THRESHOLD_FACTOR = 3.0;

%% --- 2. DATA LOADING ---
directory = uigetdir('C:\', 'Select Directory containing database.mat');
if directory == 0, return; end
load(fullfile(directory, 'database.mat'));
if istable(database), database = table2struct(database); end

results_dir = fullfile(directory, ['LBM_Analysis_3SD_' datestr(now,'yyyymmdd_HHMMSS')]);
if ~exist(results_dir, 'dir'), mkdir(results_dir); end

%% --- 3. MAIN PROCESSING LOOP ---
animalCodes = unique({database.code}, 'stable');
summaryData = {};

for ac_idx = 1:length(animalCodes)
    animalCodeStr = animalCodes{ac_idx};
    animal_indices = find(strcmp({database.code}, animalCodeStr));
    days = unique([database(animal_indices).day], 'stable');
    
    for d_idx = 1:length(days)
        curr_day = days(d_idx);
        sessionVids = animal_indices([database(animal_indices).day] == curr_day);
        
        % Data concatenation for full session
        [tempX, tempY, tempL] = deal(cell(1, 10));
        global_offset = 0; frame_indices = [];

        for i_vid = 1:length(sessionVids)
            h5_path = strtrim(database(sessionVids(i_vid)).xy);
            if ~exist(h5_path, 'file'), continue; end
            data = h5read(h5_path, '/df_with_missing/table/');
            n_frames = size(data.values_block_0, 2);
            for bp = 1:10
                cfg = BODY_PARTS_CONFIG.(bodyPartNames{bp});
                tempX{bp} = [tempX{bp}; data.values_block_0(cfg.x_row, :)'];
                tempY{bp} = [tempY{bp}; data.values_block_0(cfg.y_row, :)'];
                tempL{bp} = [tempL{bp}; data.values_block_0(cfg.L_row, :)'];
            end
            frame_indices = [frame_indices; (global_offset : global_offset + n_frames - 1)'];
            global_offset = global_offset + n_frames;
        end

        % Trimming and Interpolation
        s_trim = (SESSION_TRIM_SEC * FPS) + 1;
        e_trim = global_offset - (SESSION_TRIM_SEC * FPS);
        if s_trim > e_trim, continue; end
        trimmed_indices = frame_indices(s_trim:e_trim);
        
        speeds = NaN(length(trimmed_indices), 10);
        thresholds = NaN(1, 10);
        stats_row = {};

        for bp = 1:10
            x = tempX{bp}(s_trim:e_trim); y = tempY{bp}(s_trim:e_trim); L = tempL{bp}(s_trim:e_trim);
            x(L < LIKELIHOOD_THRESH) = NaN; y(L < LIKELIHOOD_THRESH) = NaN;
            ix = fillmissing(x, 'linear'); iy = fillmissing(y, 'linear');
            v = [NaN; sqrt(diff(ix).^2 + diff(iy).^2) * FPS];
            speeds(:, bp) = v;
            
            m_v = mean(v, 'omitnan'); s_v = std(v, 'omitnan');
            thr = m_v + (THRESHOLD_FACTOR * s_v);
            thresholds(bp) = thr;
            stats_row = [stats_row, {m_v, s_v, thr}];
        end

        % LBM DETECTION LOGIC
        % HC or TB must move; BC only counts if HC or TB are also moving
        m_hc = speeds(:, idx_hc) > thresholds(idx_hc);
        m_bc = speeds(:, idx_bc) > thresholds(idx_bc);
        m_tb = speeds(:, idx_tb) > thresholds(idx_tb);
        is_moving = (m_hc | m_tb) | (m_bc & (m_hc | m_tb));

        % Search backwards for Last Body Movement (LBM) in induction window
        lbm_sec = NaN;
        for k = length(trimmed_indices):-1:1
            abs_frame = trimmed_indices(k);
            if abs_frame >= (ISO_ADMIN_START_SEC * FPS) && abs_frame <= (LBM_WINDOW_END_SEC * FPS)
                if is_moving(k)
                    lbm_sec = abs_frame / FPS;
                    break;
                end
            end
        end

        summaryData(end+1, :) = [ {animalCodeStr, curr_day, 'T_Mean_3SD', lbm_sec}, stats_row ];
    end
end

%% --- 4. EXPORT RESULTS ---
if ~isempty(summaryData)
    headers = {'AnimalCode', 'SessionDay', 'Threshold_Type', 'Chosen_LBM_sec'};
    for b = 1:10
        headers = [headers, {['Mean_' bodyPartNames{b}], ['SD_' bodyPartNames{b}], ['Thr_' bodyPartNames{b}]}];
    end
    T = cell2table(summaryData, 'VariableNames', headers);
    writetable(T, fullfile(results_dir, 'Induction_LBM_3SD_Summary.xlsx'));
end
fprintf('Process Complete. Results saved in: %s\n', results_dir);
