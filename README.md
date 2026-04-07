# AnesthesiaDataAnalysis_MorenoGomez2026
Scripts for mouse anaesthesia emergence analysis. Supporting Moreno-Gomez M et al (2026). Includes: 1) Anaesthetized phase monitoring and 2) Emergence (FBD - First Body Movement) detection via DeepLabCut. Features consensus-based tracking (3/2 points) and kinematic filtering to identify precise awakening latency in DREADD-ligand paradigms.

=================================================================================================================================================================




==================================================================================================================================================================
# Anaesthesia Emergence Analysis Toolkit

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
