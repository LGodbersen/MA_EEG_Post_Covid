# MA_EEG_Post_Covid
Analysis of Delta and Beta EEG Power and Connectivity in Post-COVID Syndrome. This repository contains data processing, analysis scripts, and visualization tools for exploring delta and beta EEG power and connectivity in individuals experiencing post-COVID syndrome symptoms. It is a Master's Thesis based on data from a PhD project by M.Sc. Christian Neumann. This project was supervised by M.Sc. Julius Welzel.

## 1. Project Overview
Brief summary of the project, including:
Context: Describe post-COVID syndrome and its known or hypothesized neurological impact.
Objective: Explain the aim to study delta and beta EEG power and connectivity.
## 2. Data Description
- 46 participants (23 with PCS, 23 without PCS)
- age matched
- sampling rate 1000 Hz, 128 Electrodes, 3 min resting state, eyes open
- delta frequency (0.6-4 Hz) in a frontal ROI, beta frequency (14-30 Hz) in a central ROI
- high-pass (0.1 Hz) and low-pass (45 Hz) filtered, ICA decomposition and rejection with ICLabel
## 3. Analysis Methods
Signal Processing: FieldTrip Toolbox and Specparam algorithm for relative power
Connectivity Analys: imaginary part of coherence, Small World Index
Statistical Analysis: t-Tests, wilcox tests, permutation tests.
## 4. Results
Briefly summarize key findings, or provide a link to your thesis document if available.
Note any visualizations, tables, or figures used to represent findings.
## 5. Repository Structure
data/: Scripts and guidance on accessing or preprocessing EEG data.
tools/: Functions, electrode layout, MATLAB Fooof Wrapper, visualization tools.
scripts/: Main analysis scripts, organized by type of analysis (e.g., power analysis, connectivity).
results/: Code for generating visualizations, tables, and other outputs.
docs/: Additional documentation, including the full thesis.
## 6. Requirements and Setup
Toolboxes that were used:
- FieldTrip 20231127
- EEGLAB 2023.1
- BCT 2019_03_03
- if you want to fit a knee, you need the MATLAB Fooof wrapper (Python); Python needs to be installed for this
requires also the following general MATLAB toolboxes:
- Statistics and Machine Learning Toolbox
- MATLABs Optimization Toolbox
- Signal Processing Toolbox
## 7. Usage
Instructions on how to run the main analysis scripts, reproduce key findings, and generate visualizations.
## 8. Acknowledgements and References
Acknowledge any collaborators or institutions, funding sources, and relevant literature.

Data collection was conducted at the University Clinic of Schleswig-Holstein (UKSH) in Kiel at the Department of Neurology.

