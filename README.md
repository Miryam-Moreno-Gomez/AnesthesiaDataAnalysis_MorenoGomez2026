# AnesthesiaDataAnalysis_MorenoGomez2026
Scripts for mouse anaesthesia emergence analysis. Supporting Moreno-Gomez M et al (2026). Includes: 1) Anaesthetized phase monitoring and 2) Emergence (FBD - First Body Movement) detection via DeepLabCut. Features consensus-based tracking (3/2 points) and kinematic filtering to identify precise awakening latency in DREADD-ligand paradigms.

=================================================================================================================================================================

## Automated Anaesthesia Induction Analysis (LBM Detection)

### Overview
This component of the pipeline is designed to quantify the transition from a wakeful state to isoflurane-induced anaesthesia. It identifies the **Last Body Movement (LBM)**, defined as the final instance of voluntary motor activity (torso or head) before complete cessation of movement.

### Key Features
* **Dynamic Thresholding**: Implements a sensitive detection threshold based on the animal's session-specific velocity (Mean + 3.0 Standard Deviations).
* **Hierarchical Logic**: Movement is primarily anchored to the Head Centre and Tail Base. Body Centre displacement is validated only when accompanied by these anchor points to prevent tracking noise from triggering false movements.
* **Induction Window**: The analysis focuses specifically on the window following the start of isoflurane administration (from 120s to 540s).

### Output
The script generates a `Induction_LBM_3SD_Summary.xlsx` file containing:
1. **Chosen_LBM_sec**: The exact timestamp of the last voluntary movement.
2. **Body Part Statistics**: Mean speed, Standard Deviation, and the applied threshold for all 10 anatomical labels.

## Requirements
- MATLAB (R2018a+)
- Statistics and Machine Learning Toolbox
- DeepLabCut output (.h5)

### Reference
This analysis follows the methodology described in:
Moreno-Gomez et al. (2026). [cite_start]*Off-target effects of DREADD ligands revealed by an anaesthesia emergence paradigm in mice*.


==================================================================================================================================================================
# Automated Anaesthesia Emergence Analysis (FBD Detection)

This repository contains the MATLAB suite for analyzing mouse behavior during anaesthesia induction and emergence, as described in:
> **Moreno-Gomez M et al.**, *"Off-target effects of DREADD ligands revealed by an anaesthesia emergence paradigm in mice"* (2026).

## Pipeline Overview
The toolkit is divided into two specialized processing scripts to handle the different physiological states of the animal:

### Phase 1: Anaesthetized Analysis (`anaesthesia_phase_analysis.m`)
*Processes the period of deep anaesthesia.*
- Monitors stability and baseline noise floor.
- Detects subtle movements or respiratory patterns if tracked.

### Phase 2: Emergence & FBD Detection (`emergence_FBD_analysis.m`)
*Identifies the transition from anaesthesia to wakefulness.*
- **FBD (First Body Movement):** Specifically detects the first instance of coordinated movement (consensus-based) after the habituation/anaesthesia period.
- **Kinematic Filtering:** Uses displacement and path length thresholds to distinguish between "flickering" noise and true awakening latency.

## Technical Key Features
* **Hierarchical Consensus**: Validates movement only when 3 or 2 body parts (Head, Body, Tail) cross the threshold simultaneously.
* **Latency Calculation**: Automatically calculates the time from the start of the emergence trial to the FBD.
* **SVG Diagnostics**: Generates velocity profiles to visually verify the exact frame of the FBD.

## Requirements
- MATLAB (R2018a+)
- Statistics and Machine Learning Toolbox
- DeepLabCut output (.h5)

## Citation
Please cite Moreno-Gomez M et al. (2026) if you use these scripts.
