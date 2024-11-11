%% 1. Script to read in and do some general preprocessing to the data
% 0. Preliminaries
% 1. Define the current dataset
    % 1.1. Read in Markers and define trials
    % 1.2. Save the trial-definition
    % 1.3. Then define the entire data set
% 2. Read in the continuous data and apply filters (Preprocessing)
    % 2.1. High pass filter
    % 2.2. Resampling
    % 2.3 Low pass filter
    % 2.4 CAR
% 3. Safe the data for further processing

%% 0.Preliminaries
clear
close all
clc

% add source folder for functions, toolboxes, etc.
addpath(fullfile(pwd,"Src"))
% initialize fieldtrip
addpath('C:\Users\Lara Godbersen\Documents\MATLAB\fieldtrip20231127');
ft_defaults;

%% 1. Define the current dataset
proj_dir = fullfile(pwd); % automatically get path of script location, and parent dir
addpath(genpath(proj_dir)); % add dir to project to Matlab path
%indir = fullfile(proj_dir,'data\raw');% path to folder with BIDS datasets
indir = fullfile('D:\epoc_data\');% path to folder with BIDS datasets on Server/ clinic EEG: fullfile('D:\clinic_EEG\clinic EEG')
outdir = fullfile(proj_dir,'data\prep_01'); % path to prep ft data
indat = dir(indir); % content of that folder
indat = indat(startsWith({indat.name}, 'sub-')); % only keep folders that start with 'sub-' (i.e. the subjects)

%% start loop

for v = 1:length(indat) % begin to create a loop over data sets

    eegdir = fullfile(indir,indat(v).name,'eeg'); % path to the eeg folder of the current subject
    eegdat = dir(eegdir); % content of that folder
    % find file with resting state eeg data
    eegdat_eeg = eegdat(contains({eegdat.name}, 'task-restingstate_eeg.eeg')); % resting state eeg data 
    % check if resting state exists
    if isempty(eegdat_eeg);continue;end
    
    eegdat_mrk = eegdat(contains({eegdat.name}, 'restingstate_events.tsv')); % resting state eeg data events
    % extract the participant id
    tmp_id = extractBefore(eegdat_eeg.name,'_');
    
    % load the BIDS data with the preprocessing function
    cfg = [];
    cfg.dataset = fullfile(eegdat_eeg.folder,eegdat_eeg.name);
    cfg.channel = 'all';
    cfg.readbids = 'Yes';
    data_p = ft_preprocessing(cfg);
    
    % read events from events.tsv
    hdr   = ft_read_header(cfg.dataset); % get the header from the file
    event = ft_read_event(fullfile(eegdat_mrk.folder,eegdat_mrk.name), 'header', hdr, 'eventformat', 'bids_tsv'); % read in BIDS events
    
    %1.1 Read in Markers and define trial (between Experiment start and
    %Close eyes
    cfg = [];
    cfg.dataset             =  fullfile(eegdat_eeg.folder,eegdat_eeg.name);
    cfg.trialdef.eventtype  = 'Markers';
    cfg.trialfun             = 'mytrialfun';% I wrote my own trialfun -> you can find it in the Src file
    cfg = ft_definetrial(cfg);

    % 1.2. Save the trial-definition
    trl = cfg.trl;
    
    %% 2. Read in the continuous data and apply filters (Preprocessing)
    % 2.1 Highpass Filter
    % which Filter to use
    cfg.demean      = 'yes';    % remove DC offset
    cfg.hpfilter    = 'yes';
    cfg.hpfreq      = .1;       % high-pass filter, cutting everything under .1 Hz
    cfg.hpfilttype  = 'firws';
    %cfg.pad         = 'nextpow2'; cfg pad is not even a ft_preprocessing
    %option (would be cfg.padding)
    
    data_p = ft_preprocessing(cfg); 
        
    %% 2.2 Resampling 
    cfg.resamplefs = 250; % do not go below the Nyquist frequency
    cfg.method = 'resample';
    data_p = ft_resampledata(cfg, data_p);
    
    %% 2.3 Lowpass Filter    
    cfg = [];% I need to clear the configuration, otherwise it still has things like resamplefs in the cfg, which confuses ft_preprocessing
    cfg.lpfilter = 'yes'; % The same for lowpass
    cfg.lpfreq = 45; % 50Hz line noise destroys the data, only take data below
    cfg.lpfilttype = 'firws'; % Again the type
    data_clean = ft_preprocessing(cfg, data_p);

    %% 2.4 apply CAR
    % After cleaning the data, it is best to re-reference the data to
    % the average across channels to remove the influence of the
    % reference
    cfg = [];
    cfg.reref = 'yes';
    cfg.refchannel = 1:length(data_clean.label)-1; % Take all channels
    cfg.refmethod = 'avg'; % Take the average
    data_clean = ft_preprocessing(cfg,data_clean);

    %% 3. Save output
    save(fullfile(outdir,[tmp_id + "_ft_clean.mat"]),"data_clean");% save with participant ID in the name
end

