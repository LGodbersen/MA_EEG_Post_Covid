# MA_EEG_Post_Covid
Analysis of Delta and Beta EEG Power and Connectivity in Post-COVID Syndrome. This repository contains data processing, analysis scripts, and visualization tools for exploring delta and beta EEG power and connectivity in individuals experiencing post-COVID syndrome symptoms.

## 1. Project OverviewP
ost-COVID syndrome (PCS) is a condition affecting a substantial number of individuals after SARS-CoV-2 infection, often manifesting in cognitive impairments, fatigue, and other neurological symptoms (Chen et al.,2022). Existing research using electroencephalography (EEG) has highlighted irregularities in delta and beta brainwave activity in PCS patients, though the connection between these abnormalities and specific symptoms remains ambiguous (Furlanis et al.,2023; Kopanska et al.,2022,2023). In particular, previous studies have observed varied levels of delta power in PCS patients, with inconclusive findings on its relationship to cognitive symptoms (Cecchetti et al.,2022; Ortelli et al.,2023). Additionally, increased beta connectivity has been speculated to correlate with fatigue (Vecchio et al.,2017; Wu et al.,2023), a common PCS complaint.

This study aimed to examine the association between delta and beta EEG power/connectivity with self-reported cognitive difficulties and fatigue in individuals with PCS. Using high-density EEG, we investigated a range of metrics, including relative delta and beta power, the aperiodic offset and exponent, functional connectivity, and graph-based connectivity measures.
## 2. Data Description
- 46 participants (23 with PCS, 23 without PCS)
- age matched
- sampling rate 1000 Hz, 128 Electrodes, 5 min resting state, eyes open
- delta frequency (0.6-4 Hz) in a frontal ROI, beta frequency (14-30 Hz) in a central ROI
- high-pass (0.1 Hz) and low-pass (45 Hz) filtered, ICA decomposition and rejection with ICLabel
## 3. Analysis Methods
- Signal Processing: FieldTrip Toolbox and Specparam algorithm for relative power
- Connectivity Analys: imaginary part of coherence, Small World Index
- Statistical Analysis: t-Tests, wilcox tests, permutation tests.
## 4. Results
The analysis revealed no significant differences in delta or beta power, nor in delta connectivity, between PCS and control groups. However, beta functional connectivity was notably higher in the PCS group, lending support to previous research linking elevated beta connectivity with symptoms of fatigue (Buyukturkoglu et al., 2017; Vecchio et al.,2017; Wu et al.,2023). Nonetheless, further research is recommended, particularly in source space, to deepen insight into these associations and their clinical implications.
## 5. Repository Structure
- data/: tables with power and connectivity results.
- tools/: Functions, electrode layout, MATLAB Fooof Wrapper, visualization tools.
- scripts/: Main analysis scripts, organized by type of analysis (e.g., power analysis, connectivity).
- results/: Code for generating visualizations, tables, and other outputs.
- docs/: Additional documentation, including the full thesis.
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
The scripts have numbers, so this order can be used if you are working with raw data. Otherwise, the tables from the data folder can be used and loaded into the R.
## 8. Acknowledgements and References
It is a Master's Thesis based on data from a PhD project by M.Sc. Christian Neumann. This project was supervised by M.Sc. Julius Welzel and Prof. Dr. Christian Kaernbach.
Data collection was conducted at the University Clinic of Schleswig-Holstein (UKSH) in Kiel at the Department of Neurology.

Buyukturkoglu, K., Porcaro, C., Cottone, C., Cancelli, A., Inglese, M., & Tecchio, F. (2017). Simple index of functional connectivity at rest in Multiple Sclerosis fatigue. Clinical Neurophysiology, 128(5), 807–813. https://doi.org/10.1016/j.clinph.2017.02.010

Cecchetti, G., Agosta, F., Canu, E., Basaia, S., Barbieri, A., Cardamone, R., Bernasconi, M. P., Castelnovo, V., Cividini, C., Cursi, M., Vabanesi, M., Impellizzeri, M., Lazzarin, S. M., Fanelli, G. F., Minicucci, F., Giacalone, G., Falini, A., Falautano, M., Rovere-Querini, P., … Filippi, M. (2022). Cognitive, EEG, and MRI features of COVID-19 survivors: A 10-month study. Journal of Neurology, 269(7), 3400–3412. https://doi.org/10.1007/s00415-022-11047-5

Chen, C., Haupert, S. R., Zimmermann, L., Shi, X., Fritsche, L. G., & Mukherjee, B. (2022). Global Prevalence of Post-Coronavirus Disease 2019 (COVID-19) Condition or Long COVID: A Meta-Analysis and Systematic Review. The Journal of Infectious Diseases, 226(9), 1593–1607. https://doi.org/10.1093/infdis/jiac136

Furlanis, G., Buoite Stella, A., Biaduzzini, F., Bellavita, G., Frezza, N. A., Olivo, S., Menichelli, A., Lunardelli, A., Ajčević, M., & Manganotti, P. (2023). Cognitive deficit in post-acute COVID-19: An opportunity for EEG evaluation? Neurological Sciences, 44(5), 1491–1498. https://doi.org/10.1007/s10072-023-06615-0

Kopańska, M., Ochojska, D., Muchacka, R., Dejnowicz-Velitchkov, A., Banaś-Ząbczyk, A., & Szczygielski, J. (2022). Comparison of QEEG Findings before and after Onset of Post-COVID-19 Brain Fog Symptoms. Sensors, 22(17), 6606. https://doi.org/10.3390/s22176606

Kopańska, M., Rydzik, Ł., Błajda, J., Sarzyńska, I., Jachymek, K., Pałka, T., Ambroży, T., & Szczygielski, J. (2023). The Use of Quantitative Electroencephalography (QEEG) to Assess Post-COVID-19 Concentration Disorders in Professional Pilots: An Initial Concept. Brain Sciences, 13(9), 1264. https://doi.org/10.3390/brainsci13091264

Ortelli, P., Quercia, A., Cerasa, A., Dezi, S., Ferrazzoli, D., Sebastianelli, L., Saltuari, L., Versace, V., & Quartarone, A. (2023). Lowered Delta Activity in Post-COVID-19 Patients with Fatigue and Cognitive Impairment. Biomedicines, 11(8), 2228. https://doi.org/10.3390/biomedicines11082228

Vecchio, F., Miraglia, F., Porcaro, C., Cottone, C., Cancelli, A., Rossini, P. M., & Tecchio, F. (2017). Electroencephalography-Derived Sensory and Motor Network Topology in Multiple Sclerosis Fatigue. Neurorehabilitation and Neural Repair, 31(1), 56–64. https://doi.org/10.1177/1545968316656055

Wu, C.-H., De Doncker, W., & Kuppuswamy, A. (2023). Electroencephalography-Derived Functional Connectivity in Sensorimotor Networks in Post Stroke Fatigue. Brain Topography, 36(5), 727–735. https://doi.org/10.1007/s10548-023-00975-8
