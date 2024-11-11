# Overview of the script contents
## 0. Toolbox versions that I used
- FieldTrip 20231127
- EEGLAB 2023.1
- BCT 2019_03_03

## 1. read_in_and_preproc
0. Preliminaries
1. Define the current dataset
    1.1. Read in Markers and define trials
    1.2. Save the trial-definition
    1.3. Then define the entire data set
2. Read in the continuous data and apply filters (Preprocessing)
    2.1. High pass filter
    2.2. Resampling
    2.3 Low pass filter
    2.4 CAR
3. Safe the data for further processing
   
## 2. Preprocessing_main (for this you need read_in_and_preproc)
0. Preliminaries
1. read in data and convert into EEGLAB structure
2. exclude big artifacts 
3. Common average reference
4. ICA 
5. ICLabel (automatic component rejection)
6. Additional artifact removal
7. Interpolate bad channels: Repair the Cleaned and ICA-Corrected Data
8. Re-Reference 
9. Extract epochs 
10. save data in prep_power
11. Laplacian for connectivity data
12. save data in prep_connectivity
13. find out how many epochs survived
   13.1 power (5s)
   13.2 connectivity (4s)

## 3. Analysis_power
0. Preliminaries
1. Spectral Parameterization
2. create table with all relevant values
3. adding group membership and other behavioral data
4. save combined table as .csv
5. fitting knee
6. create a big table with the resulst of the python MATLAB wrapper
7. add group membership & data_behav
8. save

## 4. Visualize_power
0. Preliminaries
1. topoplots
  1.1 power
  1.2 aperiodic components
  1.3 r squared
2. plot the whole power spectrum
3. permutation test

## 5. R_Script_power
1. load packages
2. load data
3. summarise mean
4. demographics
5. outlier removal
  5.1 delta (relative and absolute)
  5.2 beta (relative)
  5.3 aperiodic components
6. export tables for topoplots
7. check requirements (normality, variances, etc.)
8. boxplots and stats
  8.1 aperiodic exponent (whole brain)
  8.2 aperiodic offset (whole brain)
  8.3 rel and abs delta frontal
  8.4 rel beta central
  8.5 tables of all EEG values
9. plot behavioral data and corr tests
  9.1 just behavioral data
  9.2 corr tests with behav - EEG data
    9.2.1 rel delta w TMTA & B-A
    9.2.2 rel delta w moca
    9.2.3 rel/abs delta w FACIT
    9.2.4 rel delta w hads
    9.2.5 rel beta w TMTA & B-A
    9.2.6 rel beta w FACIT
10. r squared
11. permutation tests

## 6. Analysis_connectivity 
0. Preliminaries
1. big loop (later ROI selection): normal, random electrodes, all electrodes
  1.1 connectivity measure delta/beta (such as coherence)
  1.2 small worldness delta/beta (threshold 0.1 bis 0.9) 
  1.3 create big table with all relevant values
  1.4 adding group membership and behavioral data
  1.5 save
2. inter frequency similarities (not pursued)

## 7. Visualize_connectivity
1. mean connectivity matrix per group
    1.1 delta
    1.2 beta
2. plots for the overview chart 
3. topoplot
4. have a look at the modularity (exploratory)

## 8. R_Script_connectivity
1. load packages
2. load data from MATLAB .csv table
3. modify dataset
  3.1 age match
  3.2 add a normalized gcc and cpl
4. have a look at the number of epochs
5. have a look at the values (SW, GCC, CPL, FC) per threshold 
6. Is the threshold itself different?
7. remove outliers
8. check requirements for icoh
9. permutation tests
10. compare normal, overall, random
11. boxplots
12. Correlations with behavioral data
  12.1 TMTA, B-A
  12.2 FACIT
13. corr of sw and coherence

