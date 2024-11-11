%% MATLAB Script for connectivity analyses of my Master's thesis
% overview
% 0. Preliminaries
% 1. big loop (later ROI selection): normal, random electrodes, all electrodes
    % 1.1 connectivity measure delta/beta (such as coherence)
    % 1.2 small worldness delta/beta (threshold 0.1 bis 0.9) 
    % 1.3 create big table with all relevant values
    % 1.4 adding group membership and behavioral data
    % 1.5 save
% 2. inter frequency similarities (not pursued)
%% 0.Preliminaries
clear
close all
clc

% initialize fieldtrip
ft_defaults;
% initialize EEGLAB
eeglab;close;
% add Brain Connectivity Toolbox
addpath(fullfile(pwd,"Src","BCT"))

%% 1. big loop (later ROI selection): normal, random electrodes, all electrodes
% set paths
proj_dir = fullfile(pwd); % automatically get path of script location, and parent dir
indir = fullfile(proj_dir,'data\prep_connectivity_01');% path to folder with BIDS datasets
outdir = fullfile(proj_dir,'data\analysis_connectivity_icoh\01'); %  alternative: 'data\analysis_connectivity_icoh'
indat = dir(indir); % content of that folder
indat = indat(startsWith({indat.name}, 'sub-')); % only keep folders that start with 'sub-' (i.e. the subjects)

% define ROIs
roi_delta = [21,104,11,39,74,38,48,81,47,19,111,24,93,92,82,91,94,95,20,49,10,58,25]; % frontal ROI (channel indices are not always the number of the actual electrode!)
roi_beta = [84,64,89,65,1,67,3,66,2,69,73,75,80,33,36,41,85,42,86,43,87,44,88,45,76,5,77,6,7,78,8,79,70,34,71,35,72]; % put channel indices here
roi_delta_random = [54,99,69,85,19,22,72,7,97,3,92,71,52,27,74,124,111,121,18,35,122,16,125];
roi_beta_random = [77,107,29,30,41,23,13,3,56,120,35,75,103,24,83,122,121,9,79,91,92,125,82,36,14,27,114,51,105,81,113,39,126,93,42,65,16];

% prepare
freq = cell(length(indat),1);
coh = cell(length(indat),1);

for m = 1:3 % three 'methods': our ROI, over whole brain, random ROI of the same size
for s = 1:length(indat)
load(fullfile(indir,indat(s).name));
data_prep = eeglab2fieldtrip(EEG_epoched_4, 'raw');% convert to FieldTrip structure

%% 1.1 connectivity measure (such as coherence)
%could look something like this: https://www.fieldtriptoolbox.org/tutorial/connectivity/
cfg           = [];
cfg.foi        = 0.3 :0.2: 30;% our high pass filter is at 0.1
cfg.method    = 'mtmfft';
cfg.taper     = 'dpss';
cfg.pad       = 'nextpow2';
cfg.output    = 'fourier';
cfg.tapsmofrq = 2;
freq{s}          = ft_freqanalysis(cfg, data_prep);

cfg           = [];
cfg.method    = 'coh';% coherence
cfg.complex   = 'imag';% compute imaginary part of coherence
coh{s}           = ft_connectivityanalysis(cfg, freq{s});

coh{s}.id =  extractBefore(indat(s).name,'_');

% connectivity matrix for delta
conn_matrix = coh{s}.cohspctrm;
conn_matrix_delta = mean(conn_matrix(:,:,3:16),3); 

% connectivity matrix for beta
conn_matrix_beta = mean(conn_matrix(:,:,58:122),3);

if m == 1
    roi_delta = roi_delta;
    roi_beta = roi_beta;
elseif m == 2
    roi_delta = roi_delta_random;
    roi_beta = roi_beta_random;
else
    roi_delta = 1:126;
    roi_beta = 1:126;
end

% I want the mean coherence per frequency (= functional connectivity measure)
conn_matrix_delta = abs(conn_matrix_delta);% making sure that the matrix is symmetrical
mean_coh_delta = mean(conn_matrix_delta(roi_delta,roi_delta,:));
coh_delta{s} = mean(mean_coh_delta);

conn_matrix_beta = abs(conn_matrix_beta); % making sure that the matrix is symmetrical
mean_coh_beta = mean(conn_matrix_beta(roi_beta,roi_beta,:));
coh_beta{s} = mean(mean_coh_beta);

%% 1.2 Small worldness
% the concept is that I have a clustering coefficient (segregation) and a
% characteristic path length (integration) for a specified ROI. Then I compute what it
% would look like if the clustering was random and devide the clustering
% that we found by the mean random clustering -> process of normalization.
% The same for characteristic path length. Then divide CC/CPL and you get a
% Small-World Index
% approach adapted from the DISCOVER EEG Pipeline
    %% 1.1 delta 
    %% threshold = 0.1
    % Threshold the connectivity matrix depending on the desired amount of edges
    sortedValues_delta_01 = sort(abs(conn_matrix_delta(:)),'descend');
    sortedValues_delta_01 = sortedValues_delta_01(~isnan(sortedValues_delta_01));
    threshold_delta_01{s} = sortedValues_delta_01(floor(0.1*length(sortedValues_delta_01)));
    
    % Binarize connectivity matrix based on the threshold (adjacency matrix)
    threshold_matrix_delta_01 = threshold_proportional(conn_matrix_delta,0.1);
    adjacency_matrix_delta_01 = weight_conversion(threshold_matrix_delta_01,'binarize');

    % extract coherence here (per threshold)
    coh_delta_01{s} = mean(threshold_matrix_delta_01(roi_delta));
    
    % maybe extract the ROI of delta here?
    adjacency_matrix_delta_roi_01 = adjacency_matrix_delta_01(roi_delta,roi_delta,:);

    % Graph analysis measures from Brain Connectivity toolbox
    % ---- Local measures ----- 
    % Degree - Number of connexions of each node
    degree_delta_01 = degrees_und(adjacency_matrix_delta_roi_01);
    % Clustering coefficient - The percentage of existing triangles surrounding
    % one node out of all possible triangles
    cc_delta_01 = clustering_coef_bu(adjacency_matrix_delta_01);
    cc_delta_01 = cc_delta_01(roi_delta);

    % ---- Global measures of segregation -----
    % Global clustering coefficient
    gcc_delta_01{s} = mean(cc_delta_01);

    % ---- Global measures of integration -----
    % Characteristic path length
    distance_delta_01 = distance_bin(adjacency_matrix_delta_01);
    distance_delta_roi_01 = distance_delta_01(roi_delta,roi_delta,:);
    [cpl_delta_01{s}, ~] = charpath(distance_delta_roi_01,0,0); % it does not include infinite distances in the calculation 
    % Global efficiency - The average of the inverse shortest path between two points of the network
    geff_delta_01{s} = efficiency_bin(adjacency_matrix_delta_roi_01); % it does include infinite distances in the calculation

    % ----- Small-worldness -----
    % Typically in small-world networks L >= Lrand but CC >> CCrand
    randN_delta_01 = makerandCIJ_und(length(adjacency_matrix_delta_01), floor(sum(sum(adjacency_matrix_delta_01)))/2);
    cc_rand_delta_01 = clustering_coef_bu(randN_delta_01);
    gcc_rand_delta_01{s} = mean(cc_rand_delta_01(roi_delta));
    distance_rand_delta_01 = distance_bin(randN_delta_01);
    distance_rand_delta_01 = distance_rand_delta_01(roi_delta,roi_delta,:);
    [cpl_rand_delta_01{s}, ~] = charpath(distance_rand_delta_01,0,0);
    smallworldness_delta_01{s} = (gcc_delta_01{s}/gcc_rand_delta_01{s}) / (cpl_delta_01{s}/cpl_rand_delta_01{s}); % I get a value here

    %% threshold = 0.2
    % Threshold the connectivity matrix depending on the desired amount of edges
    sortedValues_delta_02 = sort(abs(conn_matrix_delta(:)),'descend');
    sortedValues_delta_02 = sortedValues_delta_02(~isnan(sortedValues_delta_02));
    threshold_delta_02{s} = sortedValues_delta_02(floor(0.2*length(sortedValues_delta_02)));

    % Binarize connectivity matrix based on the threshold (adjacency matrix)
    threshold_matrix_delta_02 = threshold_proportional(conn_matrix_delta,0.2);
    adjacency_matrix_delta_02 = weight_conversion(threshold_matrix_delta_02,'binarize');
    
    % extract coherence here (per threshold)
    coh_delta_02{s} = mean(threshold_matrix_delta_02(roi_delta));
    
    % maybe extract the ROI of delta here?
    adjacency_matrix_delta_roi_02 = adjacency_matrix_delta_02(roi_delta,roi_delta,:);

    % Graph analysis measures from Brain Connectivity toolbox
    % ---- Local measures ----- 
    % Degree - Number of connexions of each node
    degree_delta_02 = degrees_und(adjacency_matrix_delta_roi_02);
    % Clustering coefficient - The percentage of existing triangles surrounding
    % one node out of all possible triangles
    cc_delta_02 = clustering_coef_bu(adjacency_matrix_delta_02);
    cc_delta_02 = cc_delta_02(roi_delta);

    % ---- Global measures of segregation -----
    % Global clustering coefficient
    gcc_delta_02{s} = mean(cc_delta_02);

    % ---- Global measures of integration -----
    % Characteristic path length
    distance_delta_02 = distance_bin(adjacency_matrix_delta_02);
    distance_delta_roi_02 = distance_delta_02(roi_delta,roi_delta,:);
    [cpl_delta_02{s}, ~] = charpath(distance_delta_roi_02,0,0); % it does not include infinite distances in the calculation 
    % Global efficiency - The average of the inverse shortest path between two points of the network
    geff_delta_02{s} = efficiency_bin(adjacency_matrix_delta_roi_02); % it does include infinite distances in the calculation

    % ----- Small-worldness -----
    % Typically in small-world networks L >= Lrand but CC >> CCrand
    randN_delta_02 = makerandCIJ_und(length(adjacency_matrix_delta_02), floor(sum(sum(adjacency_matrix_delta_02)))/2);
    cc_rand_delta_02 = clustering_coef_bu(randN_delta_02);
    gcc_rand_delta_02{s} = mean(cc_rand_delta_02(roi_delta));
    distance_rand_delta_02 = distance_bin(randN_delta_02);
    distance_rand_delta_02 = distance_rand_delta_02(roi_delta,roi_delta,:);
    [cpl_rand_delta_02{s}, ~] = charpath(distance_rand_delta_02,0,0);
    smallworldness_delta_02{s} = (gcc_delta_02{s}/gcc_rand_delta_02{s}) / (cpl_delta_02{s}/cpl_rand_delta_02{s}); % I get a value here

    %% threshold = 0.3
    % Threshold the connectivity matrix depending on the desired amount of edges
    sortedValues_delta_03 = sort(abs(conn_matrix_delta(:)),'descend');
    sortedValues_delta_03 = sortedValues_delta_03(~isnan(sortedValues_delta_03));
    threshold_delta_03{s} = sortedValues_delta_03(floor(0.3*length(sortedValues_delta_03)));

    % Binarize connectivity matrix based on the threshold (adjacency matrix)
    threshold_matrix_delta_03 = threshold_proportional(conn_matrix_delta,0.3);
    adjacency_matrix_delta_03 = weight_conversion(threshold_matrix_delta_03,'binarize');

    % extract coherence here (per threshold)
    coh_delta_03{s} = mean(threshold_matrix_delta_03(roi_delta));
    
    % maybe extract the ROI of delta here?
    adjacency_matrix_delta_roi_03 = adjacency_matrix_delta_03(roi_delta,roi_delta,:);

    % Graph analysis measures from Brain Connectivity toolbox
    % ---- Local measures ----- 
    % Degree - Number of connexions of each node
    degree_delta_03 = degrees_und(adjacency_matrix_delta_roi_03);
    % Clustering coefficient - The percentage of existing triangles surrounding
    % one node out of all possible triangles
    cc_delta_03 = clustering_coef_bu(adjacency_matrix_delta_03);
    cc_delta_03 = cc_delta_03(roi_delta);

    % ---- Global measures of segregation -----
    % Global clustering coefficient
    gcc_delta_03{s} = mean(cc_delta_03);

    % ---- Global measures of integration -----
    % Characteristic path length
    distance_delta_03 = distance_bin(adjacency_matrix_delta_03);
    distance_delta_roi_03 = distance_delta_03(roi_delta,roi_delta,:);
    [cpl_delta_03{s}, ~] = charpath(distance_delta_roi_03,0,0); % it does not include infinite distances in the calculation 
    % Global efficiency - The average of the inverse shortest path between two points of the network
    geff_delta_03{s} = efficiency_bin(adjacency_matrix_delta_roi_03); % it does include infinite distances in the calculation

    % ----- Small-worldness -----
    % Typically in small-world networks L >= Lrand but CC >> CCrand
    randN_delta_03 = makerandCIJ_und(length(adjacency_matrix_delta_03), floor(sum(sum(adjacency_matrix_delta_03)))/2);
    cc_rand_delta_03 = clustering_coef_bu(randN_delta_03);
    gcc_rand_delta_03{s} = mean(cc_rand_delta_03(roi_delta));
    distance_rand_delta_03 = distance_bin(randN_delta_03);
    distance_rand_delta_03 = distance_rand_delta_03(roi_delta,roi_delta,:);
    [cpl_rand_delta_03{s}, ~] = charpath(distance_rand_delta_03,0,0);
    smallworldness_delta_03{s} = (gcc_delta_03{s}/gcc_rand_delta_03{s}) / (cpl_delta_03{s}/cpl_rand_delta_03{s}); % I get a value here

    %% threshold = 0.4
    % Threshold the connectivity matrix depending on the desired amount of edges
    sortedValues_delta_04 = sort(abs(conn_matrix_delta(:)),'descend');
    sortedValues_delta_04 = sortedValues_delta_04(~isnan(sortedValues_delta_04));
    threshold_delta_04{s} = sortedValues_delta_04(floor(0.4*length(sortedValues_delta_04)));

    % Binarize connectivity matrix based on the threshold (adjacency matrix)
    threshold_matrix_delta_04 = threshold_proportional(conn_matrix_delta,0.4);
    adjacency_matrix_delta_04 = weight_conversion(threshold_matrix_delta_04,'binarize');

    % extract coherence here (per threshold)
    coh_delta_04{s} = mean(threshold_matrix_delta_04(roi_delta));
    
    % maybe extract the ROI of delta here?
    adjacency_matrix_delta_roi_04 = adjacency_matrix_delta_04(roi_delta,roi_delta,:);

    % Graph analysis measures from Brain Connectivity toolbox
    % ---- Local measures ----- 
    % Degree - Number of connexions of each node
    degree_delta_04 = degrees_und(adjacency_matrix_delta_roi_04);
    % Clustering coefficient - The percentage of existing triangles surrounding
    % one node out of all possible triangles
    cc_delta_04 = clustering_coef_bu(adjacency_matrix_delta_04);
    cc_delta_04 = cc_delta_04(roi_delta);

    % ---- Global measures of segregation -----
    % Global clustering coefficient
    gcc_delta_04{s} = mean(cc_delta_04);

    % ---- Global measures of integration -----
    % Characteristic path length
    distance_delta_04 = distance_bin(adjacency_matrix_delta_04);
    distance_delta_roi_04 = distance_delta_04(roi_delta,roi_delta,:);
    [cpl_delta_04{s}, ~] = charpath(distance_delta_roi_04,0,0); % it does not include infinite distances in the calculation 
    % Global efficiency - The average of the inverse shortest path between two points of the network
    geff_delta_04{s} = efficiency_bin(adjacency_matrix_delta_roi_04); % it does include infinite distances in the calculation

    % ----- Small-worldness -----
    % Typically in small-world networks L >= Lrand but CC >> CCrand
    randN_delta_04 = makerandCIJ_und(length(adjacency_matrix_delta_04), floor(sum(sum(adjacency_matrix_delta_04)))/2);
    cc_rand_delta_04 = clustering_coef_bu(randN_delta_04);
    gcc_rand_delta_04{s} = mean(cc_rand_delta_04(roi_delta));
    distance_rand_delta_04 = distance_bin(randN_delta_04);
    distance_rand_delta_04 = distance_rand_delta_04(roi_delta,roi_delta,:);
    [cpl_rand_delta_04{s}, ~] = charpath(distance_rand_delta_04,0,0);
    smallworldness_delta_04{s} = (gcc_delta_04{s}/gcc_rand_delta_04{s}) / (cpl_delta_04{s}/cpl_rand_delta_04{s}); % I get a value here

    %% threshold = 0.5
    % Threshold the connectivity matrix depending on the desired amount of edges
    sortedValues_delta_05 = sort(abs(conn_matrix_delta(:)),'descend');
    sortedValues_delta_05 = sortedValues_delta_05(~isnan(sortedValues_delta_05));
    threshold_delta_05{s} = sortedValues_delta_05(floor(0.5*length(sortedValues_delta_05)));

    % Binarize connectivity matrix based on the threshold (adjacency matrix)
    threshold_matrix_delta_05 = threshold_proportional(conn_matrix_delta,0.5);
    adjacency_matrix_delta_05 = weight_conversion(threshold_matrix_delta_05,'binarize');

    % extract coherence here (per threshold)
    coh_delta_05{s} = mean(threshold_matrix_delta_05(roi_delta));
    
    % maybe extract the ROI of delta here?
    adjacency_matrix_delta_roi_05 = adjacency_matrix_delta_05(roi_delta,roi_delta,:);

    % Graph analysis measures from Brain Connectivity toolbox
    % ---- Local measures ----- 
    % Degree - Number of connexions of each node
    degree_delta_05 = degrees_und(adjacency_matrix_delta_roi_05);
    % Clustering coefficient - The percentage of existing triangles surrounding
    % one node out of all possible triangles
    cc_delta_05 = clustering_coef_bu(adjacency_matrix_delta_05);
    cc_delta_05 = cc_delta_05(roi_delta);

    % ---- Global measures of segregation -----
    % Global clustering coefficient
    gcc_delta_05{s} = mean(cc_delta_05);

    % ---- Global measures of integration -----
    % Characteristic path length
    distance_delta_05 = distance_bin(adjacency_matrix_delta_05);
    distance_delta_roi_05 = distance_delta_05(roi_delta,roi_delta,:);
    [cpl_delta_05{s}, ~] = charpath(distance_delta_roi_05,0,0); % it does not include infinite distances in the calculation 
    % Global efficiency - The average of the inverse shortest path between two points of the network
    geff_delta_05{s} = efficiency_bin(adjacency_matrix_delta_roi_05); % it does include infinite distances in the calculation

    % ----- Small-worldness -----
    % Typically in small-world networks L >= Lrand but CC >> CCrand
    randN_delta_05 = makerandCIJ_und(length(adjacency_matrix_delta_05), floor(sum(sum(adjacency_matrix_delta_05)))/2);
    cc_rand_delta_05 = clustering_coef_bu(randN_delta_05);
    gcc_rand_delta_05{s} = mean(cc_rand_delta_05(roi_delta));
    distance_rand_delta_05 = distance_bin(randN_delta_05);
    distance_rand_delta_05 = distance_rand_delta_05(roi_delta,roi_delta,:);
    [cpl_rand_delta_05{s}, ~] = charpath(distance_rand_delta_05,0,0);
    smallworldness_delta_05{s} = (gcc_delta_05{s}/gcc_rand_delta_05{s}) / (cpl_delta_05{s}/cpl_rand_delta_05{s}); % I get a value here

    %% threshold = 0.6
    % Threshold the connectivity matrix depending on the desired amount of edges
    sortedValues_delta_06 = sort(abs(conn_matrix_delta(:)),'descend');
    sortedValues_delta_06 = sortedValues_delta_06(~isnan(sortedValues_delta_06));
    threshold_delta_06{s} = sortedValues_delta_06(floor(0.6*length(sortedValues_delta_06)));

    % Binarize connectivity matrix based on the threshold (adjacency matrix)
    threshold_matrix_delta_06 = threshold_proportional(conn_matrix_delta,0.6);
    adjacency_matrix_delta_06 = weight_conversion(threshold_matrix_delta_06,'binarize');

    % extract coherence here (per threshold)
    coh_delta_06{s} = mean(threshold_matrix_delta_06(roi_delta));
    
    % maybe extract the ROI of delta here?
    adjacency_matrix_delta_roi_06 = adjacency_matrix_delta_06(roi_delta,roi_delta,:);

    % Graph analysis measures from Brain Connectivity toolbox
    % ---- Local measures ----- 
    % Degree - Number of connexions of each node
    degree_delta_06 = degrees_und(adjacency_matrix_delta_roi_06);
    % Clustering coefficient - The percentage of existing triangles surrounding
    % one node out of all possible triangles
    cc_delta_06 = clustering_coef_bu(adjacency_matrix_delta_06);
    cc_delta_06 = cc_delta_06(roi_delta);

    % ---- Global measures of segregation -----
    % Global clustering coefficient
    gcc_delta_06{s} = mean(cc_delta_06);

    % ---- Global measures of integration -----
    % Characteristic path length
    distance_delta_06 = distance_bin(adjacency_matrix_delta_06);
    distance_delta_roi_06 = distance_delta_06(roi_delta,roi_delta,:);
    [cpl_delta_06{s}, ~] = charpath(distance_delta_roi_06,0,0); % it does not include infinite distances in the calculation 
    % Global efficiency - The average of the inverse shortest path between two points of the network
    geff_delta_06{s} = efficiency_bin(adjacency_matrix_delta_roi_06); % it does include infinite distances in the calculation

    % ----- Small-worldness -----
    % Typically in small-world networks L >= Lrand but CC >> CCrand
    randN_delta_06 = makerandCIJ_und(length(adjacency_matrix_delta_06), floor(sum(sum(adjacency_matrix_delta_06)))/2);
    cc_rand_delta_06 = clustering_coef_bu(randN_delta_06);
    gcc_rand_delta_06{s} = mean(cc_rand_delta_06(roi_delta));
    distance_rand_delta_06 = distance_bin(randN_delta_06);
    distance_rand_delta_06 = distance_rand_delta_06(roi_delta,roi_delta,:);
    [cpl_rand_delta_06{s}, ~] = charpath(distance_rand_delta_06,0,0);
    smallworldness_delta_06{s} = (gcc_delta_06{s}/gcc_rand_delta_06{s}) / (cpl_delta_06{s}/cpl_rand_delta_06{s}); % I get a value here

    %% threshold = 0.7
    % Threshold the connectivity matrix depending on the desired amount of edges
    sortedValues_delta_07 = sort(abs(conn_matrix_delta(:)),'descend');
    sortedValues_delta_07 = sortedValues_delta_07(~isnan(sortedValues_delta_07));
    threshold_delta_07{s} = sortedValues_delta_07(floor(0.7*length(sortedValues_delta_07)));

    % Binarize connectivity matrix based on the threshold (adjacency matrix)
    threshold_matrix_delta_07 = threshold_proportional(conn_matrix_delta,0.7);
    adjacency_matrix_delta_07 = weight_conversion(threshold_matrix_delta_07,'binarize');

    % extract coherence here (per threshold)
    coh_delta_07{s} = mean(threshold_matrix_delta_07(roi_delta));
    
    % maybe extract the ROI of delta here?
    adjacency_matrix_delta_roi_07 = adjacency_matrix_delta_07(roi_delta,roi_delta,:);

    % Graph analysis measures from Brain Connectivity toolbox
    % ---- Local measures ----- 
    % Degree - Number of connexions of each node
    degree_delta_07 = degrees_und(adjacency_matrix_delta_roi_07);
    % Clustering coefficient - The percentage of existing triangles surrounding
    % one node out of all possible triangles
    cc_delta_07 = clustering_coef_bu(adjacency_matrix_delta_07);
    cc_delta_07 = cc_delta_07(roi_delta);

    % ---- Global measures of segregation -----
    % Global clustering coefficient
    gcc_delta_07{s} = mean(cc_delta_07);

    % ---- Global measures of integration -----
    % Characteristic path length
    distance_delta_07 = distance_bin(adjacency_matrix_delta_07);
    distance_delta_roi_07 = distance_delta_07(roi_delta,roi_delta,:);
    [cpl_delta_07{s}, ~] = charpath(distance_delta_roi_07,0,0); % it does not include infinite distances in the calculation 
    % Global efficiency - The average of the inverse shortest path between two points of the network
    geff_delta_07{s} = efficiency_bin(adjacency_matrix_delta_roi_07); % it does include infinite distances in the calculation

    % ----- Small-worldness -----
    % Typically in small-world networks L >= Lrand but CC >> CCrand
    randN_delta_07 = makerandCIJ_und(length(adjacency_matrix_delta_07), floor(sum(sum(adjacency_matrix_delta_07)))/2);
    cc_rand_delta_07 = clustering_coef_bu(randN_delta_07);
    gcc_rand_delta_07{s} = mean(cc_rand_delta_07(roi_delta));
    distance_rand_delta_07 = distance_bin(randN_delta_07);
    distance_rand_delta_07 = distance_rand_delta_07(roi_delta,roi_delta,:);
    [cpl_rand_delta_07{s}, ~] = charpath(distance_rand_delta_07,0,0);
    smallworldness_delta_07{s} = (gcc_delta_07{s}/gcc_rand_delta_07{s}) / (cpl_delta_07{s}/cpl_rand_delta_07{s}); % I get a value here

    %% threshold = 0.8
    % Threshold the connectivity matrix depending on the desired amount of edges
    sortedValues_delta_08 = sort(abs(conn_matrix_delta(:)),'descend');
    sortedValues_delta_08 = sortedValues_delta_08(~isnan(sortedValues_delta_08));
    threshold_delta_08{s} = sortedValues_delta_08(floor(0.8*length(sortedValues_delta_08)));

    % Binarize connectivity matrix based on the threshold (adjacency matrix)
    threshold_matrix_delta_08 = threshold_proportional(conn_matrix_delta,0.8);
    adjacency_matrix_delta_08 = weight_conversion(threshold_matrix_delta_08,'binarize');

    % extract coherence here (per threshold)
    coh_delta_08{s} = mean(threshold_matrix_delta_08(roi_delta));
    
    % maybe extract the ROI of delta here?
    adjacency_matrix_delta_roi_08 = adjacency_matrix_delta_08(roi_delta,roi_delta,:);

    % Graph analysis measures from Brain Connectivity toolbox
    % ---- Local measures ----- 
    % Degree - Number of connexions of each node
    degree_delta_08 = degrees_und(adjacency_matrix_delta_roi_08);
    % Clustering coefficient - The percentage of existing triangles surrounding
    % one node out of all possible triangles
    cc_delta_08 = clustering_coef_bu(adjacency_matrix_delta_08);
    cc_delta_08 = cc_delta_08(roi_delta);

    % ---- Global measures of segregation -----
    % Global clustering coefficient
    gcc_delta_08{s} = mean(cc_delta_08);

    % ---- Global measures of integration -----
    % Characteristic path length
    distance_delta_08 = distance_bin(adjacency_matrix_delta_08);
    distance_delta_roi_08 = distance_delta_08(roi_delta,roi_delta,:);
    [cpl_delta_08{s}, ~] = charpath(distance_delta_roi_08,0,0); % it does not include infinite distances in the calculation 
    % Global efficiency - The average of the inverse shortest path between two points of the network
    geff_delta_08{s} = efficiency_bin(adjacency_matrix_delta_roi_08); % it does include infinite distances in the calculation

    % ----- Small-worldness -----
    % Typically in small-world networks L >= Lrand but CC >> CCrand
    randN_delta_08 = makerandCIJ_und(length(adjacency_matrix_delta_08), floor(sum(sum(adjacency_matrix_delta_08)))/2);
    cc_rand_delta_08 = clustering_coef_bu(randN_delta_08);
    gcc_rand_delta_08{s} = mean(cc_rand_delta_08(roi_delta));
    distance_rand_delta_08 = distance_bin(randN_delta_08);
    distance_rand_delta_08 = distance_rand_delta_08(roi_delta,roi_delta,:);
    [cpl_rand_delta_08{s}, ~] = charpath(distance_rand_delta_08,0,0);
    smallworldness_delta_08{s} = (gcc_delta_08{s}/gcc_rand_delta_08{s}) / (cpl_delta_08{s}/cpl_rand_delta_08{s}); % I get a value here

    %% threshold = 0.9
    % Threshold the connectivity matrix depending on the desired amount of edges
    sortedValues_delta_09 = sort(abs(conn_matrix_delta(:)),'descend');
    sortedValues_delta_09 = sortedValues_delta_09(~isnan(sortedValues_delta_09));
    threshold_delta_09{s} = sortedValues_delta_09(floor(0.9*length(sortedValues_delta_09)));

    % Binarize connectivity matrix based on the threshold (adjacency matrix)
    threshold_matrix_delta_09 = threshold_proportional(conn_matrix_delta,0.9);
    adjacency_matrix_delta_09 = weight_conversion(threshold_matrix_delta_09,'binarize');

    % extract coherence here (per threshold)
    coh_delta_09{s} = mean(threshold_matrix_delta_09(roi_delta));
    
    % maybe extract the ROI of delta here?
    adjacency_matrix_delta_roi_09 = adjacency_matrix_delta_09(roi_delta,roi_delta,:);

    % Graph analysis measures from Brain Connectivity toolbox
    % ---- Local measures ----- 
    % Degree - Number of connexions of each node
    degree_delta_09 = degrees_und(adjacency_matrix_delta_roi_09);
    % Clustering coefficient - The percentage of existing triangles surrounding
    % one node out of all possible triangles
    cc_delta_09 = clustering_coef_bu(adjacency_matrix_delta_09);
    cc_delta_09 = cc_delta_09(roi_delta);

    % ---- Global measures of segregation -----
    % Global clustering coefficient
    gcc_delta_09{s} = mean(cc_delta_09);

    % ---- Global measures of integration -----
    % Characteristic path length
    distance_delta_09 = distance_bin(adjacency_matrix_delta_09);
    distance_delta_roi_09 = distance_delta_09(roi_delta,roi_delta,:);
    [cpl_delta_09{s}, ~] = charpath(distance_delta_roi_09,0,0); % it does not include infinite distances in the calculation 
    % Global efficiency - The average of the inverse shortest path between two points of the network
    geff_delta_09{s} = efficiency_bin(adjacency_matrix_delta_roi_09); % it does include infinite distances in the calculation

    % ----- Small-worldness -----
    % Typically in small-world networks L >= Lrand but CC >> CCrand
    randN_delta_09 = makerandCIJ_und(length(adjacency_matrix_delta_09), floor(sum(sum(adjacency_matrix_delta_09)))/2);
    cc_rand_delta_09 = clustering_coef_bu(randN_delta_09);
    gcc_rand_delta_09{s} = mean(cc_rand_delta_09(roi_delta));
    distance_rand_delta_09 = distance_bin(randN_delta_09);
    distance_rand_delta_09 = distance_rand_delta_09(roi_delta,roi_delta,:);
    [cpl_rand_delta_09{s}, ~] = charpath(distance_rand_delta_09,0,0);
    smallworldness_delta_09{s} = (gcc_delta_09{s}/gcc_rand_delta_09{s}) / (cpl_delta_09{s}/cpl_rand_delta_09{s}); % I get a value here

    %% let's do beta
    %% threshold = 0.1
    % Threshold the connectivity matrix depending on the desired amount of edges
    sortedValues_beta_01 = sort(abs(conn_matrix_beta(:)),'descend');
    sortedValues_beta_01 = sortedValues_beta_01(~isnan(sortedValues_beta_01));
    threshold_beta_01{s} = sortedValues_beta_01(floor(0.1*length(sortedValues_beta_01)));

    % Binarize connectivity matrix based on the threshold (adjacency matrix)
    threshold_matrix_beta_01 = threshold_proportional(conn_matrix_beta,0.1);
    adjacency_matrix_beta_01 = weight_conversion(threshold_matrix_beta_01,'binarize');

    % extract coherence here (per threshold)
    coh_beta_01{s} = mean(threshold_matrix_beta_01(roi_beta));
    
    % maybe extract the ROI of beta here?
    adjacency_matrix_beta_roi_01 = adjacency_matrix_beta_01(roi_beta,roi_beta,:);

    % Graph analysis measures from Brain Connectivity toolbox
    % ---- Local measures ----- 
    % Degree - Number of connexions of each node
    degree_beta_01 = degrees_und(adjacency_matrix_beta_roi_01);
    % Clustering coefficient - The percentage of existing triangles surrounding
    % one node out of all possible triangles
    cc_beta_01 = clustering_coef_bu(adjacency_matrix_beta_01);
    cc_beta_01 = cc_beta_01(roi_beta);

    % ---- Global measures of segregation -----
    % Global clustering coefficient
    gcc_beta_01{s} = mean(cc_beta_01);

    % ---- Global measures of integration -----
    % Characteristic path length
    distance_beta_01 = distance_bin(adjacency_matrix_beta_01);
    distance_beta_roi_01 = distance_beta_01(roi_beta,roi_beta,:);
    [cpl_beta_01{s}, ~] = charpath(distance_beta_roi_01,0,0); % it does not include infinite distances in the calculation 
    % Global efficiency - The average of the inverse shortest path between two points of the network
    geff_beta_01{s} = efficiency_bin(adjacency_matrix_beta_roi_01); % it does include infinite distances in the calculation

    % ----- Small-worldness -----
    % Typically in small-world networks L >= Lrand but CC >> CCrand
    randN_beta_01 = makerandCIJ_und(length(adjacency_matrix_beta_01), floor(sum(sum(adjacency_matrix_beta_01)))/2);
    cc_rand_beta_01 = clustering_coef_bu(randN_beta_01);
    gcc_rand_beta_01{s} = mean(cc_rand_beta_01(roi_beta));
    distance_rand_beta_01 = distance_bin(randN_beta_01);
    distance_rand_beta_01 = distance_rand_beta_01(roi_beta,roi_beta,:);
    [cpl_rand_beta_01{s}, ~] = charpath(distance_rand_beta_01,0,0);
    smallworldness_beta_01{s} = (gcc_beta_01{s}/gcc_rand_beta_01{s}) / (cpl_beta_01{s}/cpl_rand_beta_01{s}); % I get a value here

    %% threshold = 0.2
    % Threshold the connectivity matrix depending on the desired amount of edges
    sortedValues_beta_02 = sort(abs(conn_matrix_beta(:)),'descend');
    sortedValues_beta_02 = sortedValues_beta_02(~isnan(sortedValues_beta_02));
    threshold_beta_02{s} = sortedValues_beta_02(floor(0.2*length(sortedValues_beta_02)));

    % Binarize connectivity matrix based on the threshold (adjacency matrix)
    threshold_matrix_beta_02 = threshold_proportional(conn_matrix_beta,0.2);
    adjacency_matrix_beta_02 = weight_conversion(threshold_matrix_beta_02,'binarize');

    % extract coherence here (per threshold)
    coh_beta_02{s} = mean(threshold_matrix_beta_02(roi_beta));
    
    % maybe extract the ROI of beta here?
    adjacency_matrix_beta_roi_02 = adjacency_matrix_beta_02(roi_beta,roi_beta,:);

    % Graph analysis measures from Brain Connectivity toolbox
    % ---- Local measures ----- 
    % Degree - Number of connexions of each node
    degree_beta_02 = degrees_und(adjacency_matrix_beta_roi_02);
    % Clustering coefficient - The percentage of existing triangles surrounding
    % one node out of all possible triangles
    cc_beta_02 = clustering_coef_bu(adjacency_matrix_beta_02);
    cc_beta_02 = cc_beta_02(roi_beta);

    % ---- Global measures of segregation -----
    % Global clustering coefficient
    gcc_beta_02{s} = mean(cc_beta_02);

    % ---- Global measures of integration -----
    % Characteristic path length
    distance_beta_02 = distance_bin(adjacency_matrix_beta_02);
    distance_beta_roi_02 = distance_beta_02(roi_beta,roi_beta,:);
    [cpl_beta_02{s}, ~] = charpath(distance_beta_roi_02,0,0); % it does not include infinite distances in the calculation 
    % Global efficiency - The average of the inverse shortest path between two points of the network
    geff_beta_02{s} = efficiency_bin(adjacency_matrix_beta_roi_02); % it does include infinite distances in the calculation

    % ----- Small-worldness -----
    % Typically in small-world networks L >= Lrand but CC >> CCrand
    randN_beta_02 = makerandCIJ_und(length(adjacency_matrix_beta_02), floor(sum(sum(adjacency_matrix_beta_02)))/2);
    cc_rand_beta_02 = clustering_coef_bu(randN_beta_02);
    gcc_rand_beta_02{s} = mean(cc_rand_beta_02(roi_beta));
    distance_rand_beta_02 = distance_bin(randN_beta_02);
    distance_rand_beta_02 = distance_rand_beta_02(roi_beta,roi_beta,:);
    [cpl_rand_beta_02{s}, ~] = charpath(distance_rand_beta_02,0,0);
    smallworldness_beta_02{s} = (gcc_beta_02{s}/gcc_rand_beta_02{s}) / (cpl_beta_02{s}/cpl_rand_beta_02{s}); % I get a value here

    %% threshold = 0.3
    % Threshold the connectivity matrix depending on the desired amount of edges
    sortedValues_beta_03 = sort(abs(conn_matrix_beta(:)),'descend');
    sortedValues_beta_03 = sortedValues_beta_03(~isnan(sortedValues_beta_03));
    threshold_beta_03{s} = sortedValues_beta_03(floor(0.3*length(sortedValues_beta_03)));

    % Binarize connectivity matrix based on the threshold (adjacency matrix)
    threshold_matrix_beta_03 = threshold_proportional(conn_matrix_beta,0.3);
    adjacency_matrix_beta_03 = weight_conversion(threshold_matrix_beta_03,'binarize');

    % extract coherence here (per threshold)
    coh_beta_03{s} = mean(threshold_matrix_beta_03(roi_beta));
    
    % maybe extract the ROI of beta here?
    adjacency_matrix_beta_roi_03 = adjacency_matrix_beta_03(roi_beta,roi_beta,:);

    % Graph analysis measures from Brain Connectivity toolbox
    % ---- Local measures ----- 
    % Degree - Number of connexions of each node
    degree_beta_03 = degrees_und(adjacency_matrix_beta_roi_03);
    % Clustering coefficient - The percentage of existing triangles surrounding
    % one node out of all possible triangles
    cc_beta_03 = clustering_coef_bu(adjacency_matrix_beta_03);
    cc_beta_03 = cc_beta_03(roi_beta);

    % ---- Global measures of segregation -----
    % Global clustering coefficient
    gcc_beta_03{s} = mean(cc_beta_03);

    % ---- Global measures of integration -----
    % Characteristic path length
    distance_beta_03 = distance_bin(adjacency_matrix_beta_03);
    distance_beta_roi_03 = distance_beta_03(roi_beta,roi_beta,:);
    [cpl_beta_03{s}, ~] = charpath(distance_beta_roi_03,0,0); % it does not include infinite distances in the calculation 
    % Global efficiency - The average of the inverse shortest path between two points of the network
    geff_beta_03{s} = efficiency_bin(adjacency_matrix_beta_roi_03); % it does include infinite distances in the calculation

    % ----- Small-worldness -----
    % Typically in small-world networks L >= Lrand but CC >> CCrand
    randN_beta_03 = makerandCIJ_und(length(adjacency_matrix_beta_03), floor(sum(sum(adjacency_matrix_beta_03)))/2);
    cc_rand_beta_03 = clustering_coef_bu(randN_beta_03);
    gcc_rand_beta_03{s} = mean(cc_rand_beta_03(roi_beta));
    distance_rand_beta_03 = distance_bin(randN_beta_03);
    distance_rand_beta_03 = distance_rand_beta_03(roi_beta,roi_beta,:);
    [cpl_rand_beta_03{s}, ~] = charpath(distance_rand_beta_03,0,0);
    smallworldness_beta_03{s} = (gcc_beta_03{s}/gcc_rand_beta_03{s}) / (cpl_beta_03{s}/cpl_rand_beta_03{s}); % I get a value here

    %% threshold = 0.4
    % Threshold the connectivity matrix depending on the desired amount of edges
    sortedValues_beta_04 = sort(abs(conn_matrix_beta(:)),'descend');
    sortedValues_beta_04 = sortedValues_beta_04(~isnan(sortedValues_beta_04));
    threshold_beta_04{s} = sortedValues_beta_04(floor(0.4*length(sortedValues_beta_04)));

    % Binarize connectivity matrix based on the threshold (adjacency matrix)
    threshold_matrix_beta_04 = threshold_proportional(conn_matrix_beta,0.4);
    adjacency_matrix_beta_04 = weight_conversion(threshold_matrix_beta_04,'binarize');

    % extract coherence here (per threshold)
    coh_beta_04{s} = mean(threshold_matrix_beta_04(roi_beta));
    
    % maybe extract the ROI of beta here?
    adjacency_matrix_beta_roi_04 = adjacency_matrix_beta_04(roi_beta,roi_beta,:);

    % Graph analysis measures from Brain Connectivity toolbox
    % ---- Local measures ----- 
    % Degree - Number of connexions of each node
    degree_beta_04 = degrees_und(adjacency_matrix_beta_roi_04);
    % Clustering coefficient - The percentage of existing triangles surrounding
    % one node out of all possible triangles
    cc_beta_04 = clustering_coef_bu(adjacency_matrix_beta_04);
    cc_beta_04 = cc_beta_04(roi_beta);

    % ---- Global measures of segregation -----
    % Global clustering coefficient
    gcc_beta_04{s} = mean(cc_beta_04);

    % ---- Global measures of integration -----
    % Characteristic path length
    distance_beta_04 = distance_bin(adjacency_matrix_beta_04);
    distance_beta_roi_04 = distance_beta_04(roi_beta,roi_beta,:);
    [cpl_beta_04{s}, ~] = charpath(distance_beta_roi_04,0,0); % it does not include infinite distances in the calculation 
    % Global efficiency - The average of the inverse shortest path between two points of the network
    geff_beta_04{s} = efficiency_bin(adjacency_matrix_beta_roi_04); % it does include infinite distances in the calculation

    % ----- Small-worldness -----
    % Typically in small-world networks L >= Lrand but CC >> CCrand
    randN_beta_04 = makerandCIJ_und(length(adjacency_matrix_beta_04), floor(sum(sum(adjacency_matrix_beta_04)))/2);
    cc_rand_beta_04 = clustering_coef_bu(randN_beta_04);
    gcc_rand_beta_04{s} = mean(cc_rand_beta_04(roi_beta));
    distance_rand_beta_04 = distance_bin(randN_beta_04);
    distance_rand_beta_04 = distance_rand_beta_04(roi_beta,roi_beta,:);
    [cpl_rand_beta_04{s}, ~] = charpath(distance_rand_beta_04,0,0);
    smallworldness_beta_04{s} = (gcc_beta_04{s}/gcc_rand_beta_04{s}) / (cpl_beta_04{s}/cpl_rand_beta_04{s}); % I get a value here

    %% threshold = 0.5
    % Threshold the connectivity matrix depending on the desired amount of edges
    sortedValues_beta_05 = sort(abs(conn_matrix_beta(:)),'descend');
    sortedValues_beta_05 = sortedValues_beta_05(~isnan(sortedValues_beta_05));
    threshold_beta_05{s} = sortedValues_beta_05(floor(0.5*length(sortedValues_beta_05)));

    % Binarize connectivity matrix based on the threshold (adjacency matrix)
    threshold_matrix_beta_05 = threshold_proportional(conn_matrix_beta,0.5);
    adjacency_matrix_beta_05 = weight_conversion(threshold_matrix_beta_05,'binarize');

    % extract coherence here (per threshold)
    coh_beta_05{s} = mean(threshold_matrix_beta_05(roi_beta));
    
    % maybe extract the ROI of beta here?
    adjacency_matrix_beta_roi_05 = adjacency_matrix_beta_05(roi_beta,roi_beta,:);

    % Graph analysis measures from Brain Connectivity toolbox
    % ---- Local measures ----- 
    % Degree - Number of connexions of each node
    degree_beta_05 = degrees_und(adjacency_matrix_beta_roi_05);
    % Clustering coefficient - The percentage of existing triangles surrounding
    % one node out of all possible triangles
    cc_beta_05 = clustering_coef_bu(adjacency_matrix_beta_05);
    cc_beta_05 = cc_beta_05(roi_beta);

    % ---- Global measures of segregation -----
    % Global clustering coefficient
    gcc_beta_05{s} = mean(cc_beta_05);

    % ---- Global measures of integration -----
    % Characteristic path length
    distance_beta_05 = distance_bin(adjacency_matrix_beta_05);
    distance_beta_roi_05 = distance_beta_05(roi_beta,roi_beta,:);
    [cpl_beta_05{s}, ~] = charpath(distance_beta_roi_05,0,0); % it does not include infinite distances in the calculation 
    % Global efficiency - The average of the inverse shortest path between two points of the network
    geff_beta_05{s} = efficiency_bin(adjacency_matrix_beta_roi_05); % it does include infinite distances in the calculation

    % ----- Small-worldness -----
    % Typically in small-world networks L >= Lrand but CC >> CCrand
    randN_beta_05 = makerandCIJ_und(length(adjacency_matrix_beta_05), floor(sum(sum(adjacency_matrix_beta_05)))/2);
    cc_rand_beta_05 = clustering_coef_bu(randN_beta_05);
    gcc_rand_beta_05{s} = mean(cc_rand_beta_05(roi_beta));
    distance_rand_beta_05 = distance_bin(randN_beta_05);
    distance_rand_beta_05 = distance_rand_beta_05(roi_beta,roi_beta,:);
    [cpl_rand_beta_05{s}, ~] = charpath(distance_rand_beta_05,0,0);
    smallworldness_beta_05{s} = (gcc_beta_05{s}/gcc_rand_beta_05{s}) / (cpl_beta_05{s}/cpl_rand_beta_05{s}); % I get a value here

    %% threshold = 0.6
    % Threshold the connectivity matrix depending on the desired amount of edges
    sortedValues_beta_06 = sort(abs(conn_matrix_beta(:)),'descend');
    sortedValues_beta_06 = sortedValues_beta_06(~isnan(sortedValues_beta_06));
    threshold_beta_06{s} = sortedValues_beta_06(floor(0.6*length(sortedValues_beta_06)));

    % Binarize connectivity matrix based on the threshold (adjacency matrix)
    threshold_matrix_beta_06 = threshold_proportional(conn_matrix_beta,0.6);
    adjacency_matrix_beta_06 = weight_conversion(threshold_matrix_beta_06,'binarize');

    % extract coherence here (per threshold)
    coh_beta_06{s} = mean(threshold_matrix_beta_06(roi_beta));
    
    % maybe extract the ROI of beta here?
    adjacency_matrix_beta_roi_06 = adjacency_matrix_beta_06(roi_beta,roi_beta,:);

    % Graph analysis measures from Brain Connectivity toolbox
    % ---- Local measures ----- 
    % Degree - Number of connexions of each node
    degree_beta_06 = degrees_und(adjacency_matrix_beta_roi_06);
    % Clustering coefficient - The percentage of existing triangles surrounding
    % one node out of all possible triangles
    cc_beta_06 = clustering_coef_bu(adjacency_matrix_beta_06);
    cc_beta_06 = cc_beta_06(roi_beta);

    % ---- Global measures of segregation -----
    % Global clustering coefficient
    gcc_beta_06{s} = mean(cc_beta_06);

    % ---- Global measures of integration -----
    % Characteristic path length
    distance_beta_06 = distance_bin(adjacency_matrix_beta_06);
    distance_beta_roi_06 = distance_beta_06(roi_beta,roi_beta,:);
    [cpl_beta_06{s}, ~] = charpath(distance_beta_roi_06,0,0); % it does not include infinite distances in the calculation 
    % Global efficiency - The average of the inverse shortest path between two points of the network
    geff_beta_06{s} = efficiency_bin(adjacency_matrix_beta_roi_06); % it does include infinite distances in the calculation

    % ----- Small-worldness -----
    % Typically in small-world networks L >= Lrand but CC >> CCrand
    randN_beta_06 = makerandCIJ_und(length(adjacency_matrix_beta_06), floor(sum(sum(adjacency_matrix_beta_06)))/2);
    cc_rand_beta_06 = clustering_coef_bu(randN_beta_06);
    gcc_rand_beta_06{s} = mean(cc_rand_beta_06(roi_beta));
    distance_rand_beta_06 = distance_bin(randN_beta_06);
    distance_rand_beta_06 = distance_rand_beta_06(roi_beta,roi_beta,:);
    [cpl_rand_beta_06{s}, ~] = charpath(distance_rand_beta_06,0,0);
    smallworldness_beta_06{s} = (gcc_beta_06{s}/gcc_rand_beta_06{s}) / (cpl_beta_06{s}/cpl_rand_beta_06{s}); % I get a value here

    %% threshold = 0.7
    % Threshold the connectivity matrix depending on the desired amount of edges
    sortedValues_beta_07 = sort(abs(conn_matrix_beta(:)),'descend');
    sortedValues_beta_07 = sortedValues_beta_07(~isnan(sortedValues_beta_07));
    threshold_beta_07{s} = sortedValues_beta_07(floor(0.7*length(sortedValues_beta_07)));

    % Binarize connectivity matrix based on the threshold (adjacency matrix)
    threshold_matrix_beta_07 = threshold_proportional(conn_matrix_beta,0.7);
    adjacency_matrix_beta_07 = weight_conversion(threshold_matrix_beta_07,'binarize');

    % extract coherence here (per threshold)
    coh_beta_07{s} = mean(threshold_matrix_beta_07(roi_beta));
    
    % maybe extract the ROI of beta here?
    adjacency_matrix_beta_roi_07 = adjacency_matrix_beta_07(roi_beta,roi_beta,:);

    % Graph analysis measures from Brain Connectivity toolbox
    % ---- Local measures ----- 
    % Degree - Number of connexions of each node
    degree_beta_07 = degrees_und(adjacency_matrix_beta_roi_07);
    % Clustering coefficient - The percentage of existing triangles surrounding
    % one node out of all possible triangles
    cc_beta_07 = clustering_coef_bu(adjacency_matrix_beta_07);
    cc_beta_07 = cc_beta_07(roi_beta);

    % ---- Global measures of segregation -----
    % Global clustering coefficient
    gcc_beta_07{s} = mean(cc_beta_07);

    % ---- Global measures of integration -----
    % Characteristic path length
    distance_beta_07 = distance_bin(adjacency_matrix_beta_07);
    distance_beta_roi_07 = distance_beta_07(roi_beta,roi_beta,:);
    [cpl_beta_07{s}, ~] = charpath(distance_beta_roi_07,0,0); % it does not include infinite distances in the calculation 
    % Global efficiency - The average of the inverse shortest path between two points of the network
    geff_beta_07{s} = efficiency_bin(adjacency_matrix_beta_roi_07); % it does include infinite distances in the calculation

    % ----- Small-worldness -----
    % Typically in small-world networks L >= Lrand but CC >> CCrand
    randN_beta_07 = makerandCIJ_und(length(adjacency_matrix_beta_07), floor(sum(sum(adjacency_matrix_beta_07)))/2);
    cc_rand_beta_07 = clustering_coef_bu(randN_beta_07);
    gcc_rand_beta_07{s} = mean(cc_rand_beta_07(roi_beta));
    distance_rand_beta_07 = distance_bin(randN_beta_07);
    distance_rand_beta_07 = distance_rand_beta_07(roi_beta,roi_beta,:);
    [cpl_rand_beta_07{s}, ~] = charpath(distance_rand_beta_07,0,0);
    smallworldness_beta_07{s} = (gcc_beta_07{s}/gcc_rand_beta_07{s}) / (cpl_beta_07{s}/cpl_rand_beta_07{s}); % I get a value here

    %% threshold = 0.8
    % Threshold the connectivity matrix depending on the desired amount of edges
    sortedValues_beta_08 = sort(abs(conn_matrix_beta(:)),'descend');
    sortedValues_beta_08 = sortedValues_beta_08(~isnan(sortedValues_beta_08));
    threshold_beta_08{s} = sortedValues_beta_08(floor(0.8*length(sortedValues_beta_08)));

    % Binarize connectivity matrix based on the threshold (adjacency matrix)
    threshold_matrix_beta_08 = threshold_proportional(conn_matrix_beta,0.8);
    adjacency_matrix_beta_08 = weight_conversion(threshold_matrix_beta_08,'binarize');

    % extract coherence here (per threshold)
    coh_beta_08{s} = mean(threshold_matrix_beta_08(roi_beta));
    
    % maybe extract the ROI of beta here?
    adjacency_matrix_beta_roi_08 = adjacency_matrix_beta_08(roi_beta,roi_beta,:);

    % Graph analysis measures from Brain Connectivity toolbox
    % ---- Local measures ----- 
    % Degree - Number of connexions of each node
    degree_beta_08 = degrees_und(adjacency_matrix_beta_roi_08);
    % Clustering coefficient - The percentage of existing triangles surrounding
    % one node out of all possible triangles
    cc_beta_08 = clustering_coef_bu(adjacency_matrix_beta_08);
    cc_beta_08 = cc_beta_08(roi_beta);

    % ---- Global measures of segregation -----
    % Global clustering coefficient
    gcc_beta_08{s} = mean(cc_beta_08);

    % ---- Global measures of integration -----
    % Characteristic path length
    distance_beta_08 = distance_bin(adjacency_matrix_beta_08);
    distance_beta_roi_08 = distance_beta_08(roi_beta,roi_beta,:);
    [cpl_beta_08{s}, ~] = charpath(distance_beta_roi_08,0,0); % it does not include infinite distances in the calculation 
    % Global efficiency - The average of the inverse shortest path between two points of the network
    geff_beta_08{s} = efficiency_bin(adjacency_matrix_beta_roi_08); % it does include infinite distances in the calculation

    % ----- Small-worldness -----
    % Typically in small-world networks L >= Lrand but CC >> CCrand
    randN_beta_08 = makerandCIJ_und(length(adjacency_matrix_beta_08), floor(sum(sum(adjacency_matrix_beta_08)))/2);
    cc_rand_beta_08 = clustering_coef_bu(randN_beta_08);
    gcc_rand_beta_08{s} = mean(cc_rand_beta_08(roi_beta));
    distance_rand_beta_08 = distance_bin(randN_beta_08);
    distance_rand_beta_08 = distance_rand_beta_08(roi_beta,roi_beta,:);
    [cpl_rand_beta_08{s}, ~] = charpath(distance_rand_beta_08,0,0);
    smallworldness_beta_08{s} = (gcc_beta_08{s}/gcc_rand_beta_08{s}) / (cpl_beta_08{s}/cpl_rand_beta_08{s}); % I get a value here

    %% threshold = 0.9
    % Threshold the connectivity matrix depending on the desired amount of edges
    sortedValues_beta_09 = sort(abs(conn_matrix_beta(:)),'descend');
    sortedValues_beta_09 = sortedValues_beta_09(~isnan(sortedValues_beta_09));
    threshold_beta_09{s} = sortedValues_beta_09(floor(0.9*length(sortedValues_beta_09)));

    % Binarize connectivity matrix based on the threshold (adjacency matrix)
    threshold_matrix_beta_09 = threshold_proportional(conn_matrix_beta,0.9);
    adjacency_matrix_beta_09 = weight_conversion(threshold_matrix_beta_09,'binarize');

    % extract coherence here (per threshold)
    coh_beta_09{s} = mean(threshold_matrix_beta_09(roi_beta));
    
    % maybe extract the ROI of beta here?
    adjacency_matrix_beta_roi_09 = adjacency_matrix_beta_09(roi_beta,roi_beta,:);

    % Graph analysis measures from Brain Connectivity toolbox
    % ---- Local measures ----- 
    % Degree - Number of connexions of each node
    degree_beta_09 = degrees_und(adjacency_matrix_beta_roi_09);
    % Clustering coefficient - The percentage of existing triangles surrounding
    % one node out of all possible triangles
    cc_beta_09 = clustering_coef_bu(adjacency_matrix_beta_09);
    cc_beta_09 = cc_beta_09(roi_beta);

    % ---- Global measures of segregation -----
    % Global clustering coefficient
    gcc_beta_09{s} = mean(cc_beta_09);

    % ---- Global measures of integration -----
    % Characteristic path length
    distance_beta_09 = distance_bin(adjacency_matrix_beta_09);
    distance_beta_roi_09 = distance_beta_09(roi_beta,roi_beta,:);
    [cpl_beta_09{s}, ~] = charpath(distance_beta_roi_09,0,0); % it does not include infinite distances in the calculation 
    % Global efficiency - The average of the inverse shortest path between two points of the network
    geff_beta_09{s} = efficiency_bin(adjacency_matrix_beta_roi_09); % it does include infinite distances in the calculation

    % ----- Small-worldness -----
    % Typically in small-world networks L >= Lrand but CC >> CCrand
    randN_beta_09 = makerandCIJ_und(length(adjacency_matrix_beta_09), floor(sum(sum(adjacency_matrix_beta_09)))/2);
    cc_rand_beta_09 = clustering_coef_bu(randN_beta_09);
    gcc_rand_beta_09{s} = mean(cc_rand_beta_09(roi_beta));
    distance_rand_beta_09 = distance_bin(randN_beta_09);
    distance_rand_beta_09 = distance_rand_beta_09(roi_beta,roi_beta,:);
    [cpl_rand_beta_09{s}, ~] = charpath(distance_rand_beta_09,0,0);
    smallworldness_beta_09{s} = (gcc_beta_09{s}/gcc_rand_beta_09{s}) / (cpl_beta_09{s}/cpl_rand_beta_09{s}); % I get a value here

end

%% 1.3 Create a big table with all relevant values (different ROI extraction)
% Initialize an empty table
bigTable = table();
run(fullfile(pwd,'\Src\generate_big_table.m'));  % Execute generate_big_table.m is in Src folder

   
   %% 1.4 adding group membership and other behavioral data
data_behav = readtable(fullfile(proj_dir, "\data\PuG\participants.tsv"), "FileType","text",'Delimiter', '\t');
bigTable_combined = outerjoin(bigTable, data_behav,'MergeKeys',true);

   %% 1.5 save
if m == 1
    csvFile = 'table_graph_measures_ROI.csv';
    writetable(bigTable_combined, fullfile(outdir,csvFile));
elseif m == 2
    csvFile = 'table_graph_measures_ROI_random.csv';
    writetable(bigTable_combined, fullfile(outdir,csvFile));
else
    csvFile = 'table_graph_measures_ROI_overall.csv';
    writetable(bigTable_combined, fullfile(outdir,csvFile));
end

end


%% 2. inter-frequency (not pursued further)
% for this there is a github page with functions to compute the so-called
% MPC, multi-participation coefficient (https://github.com/brain-network/bnt/tree/develop)

% add BNT Toolbox
%addpath('C:\Users\Lara Godbersen\Documents\MATLAB\BNT_Toolbox\bnt-develop');

% but in order to make this supramatrix thing I need a tensor (channel x
% time x frequency)