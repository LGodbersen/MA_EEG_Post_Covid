%% Script to visualize the power data (topoplots and power spectrum)
% overview
% 0. Preliminaries
% 1. topoplots
    % 1.1 power
    % 1.2 aperiodic components
    % 1.3 r squared
% 2. plot the whole power spectrum
% 3. permutation test
%% 0. Preliminaries

clear all
close all
clc

proj_dir = fullfile(pwd);
indir = fullfile(proj_dir,'data\analysis_power');
eeglab;close;% initialize EEGLAB, you will need it for the plots


roi_delta = [21,102,11,37,72,36,46,79,45,19,109,24,91,90,80,89,92,93,20,47,10,56,25]; % frontal ROI (channel indices are not always the number of the actual electrode!)
roi_beta = [82,31,62,34,87,63,1,65,3,64,2,67,71,73,78,31,34,39,83,40,84,41,85,42,86,43,74,5,75,6,7,76,8,77,68,32,69,33,70]; % put channel indices here

%% 1. topoplots
%% 1.1 power
% take a data vector and plot it
% for this to work you need to load in one dataset in order to have the
% structure of the layout

data_beta_pcs = importdata('data\analysis_power\export_beta_pcs.txt');
data_beta_c = importdata('data\analysis_power\export_beta_c.txt');
data_delta_pcs = importdata('data\analysis_power\export_delta_pcs.txt');
data_delta_c = importdata('data\analysis_power\export_delta_c.txt');

% you need to load one preprocessed data set for this in order to get
% EEG_epoched_5.chanlocs

  f =  figure;
topoplot(data_beta_pcs.data,EEG_epoched_5.chanlocs,'colormap',viridis,'electrodes','on','emarker2', {roi_beta,'o','w',3,1})
hc=colorbar;
caxis([0.346 2.97])
xlabel(hc,'beta Power [μV^2]');
%title ('mit PCS ');
set(findobj(gca,'type','patch'),'facecolor', '#F59541'); % Change [0.5, 0.5, 0.5] to your desired color
set(gca, 'FontSize', 17);
save_fig(f,'C:\Users\Lara Godbersen\Documents\GitHub\Masters-thesis\Plots\Final\','beta_pcs', 'fontsize',17, 'figsize', [0 0 10 6.5]);


 f = figure;
topoplot(data_beta_c.data,EEG_epoched_5.chanlocs,'colormap',viridis,'electrodes','on','emarker2', {roi_beta,'o','w',3,1})
hc=colorbar;
caxis([0.346 2.97])
xlabel(hc,'beta Power [μV^2]');
%title ('ohne PCS ');
set(findobj(gca,'type','patch'),'facecolor', '#02CAF5'); % Change [0.5, 0.5, 0.5] to your desired color
set(gca, 'FontSize', 17);
save_fig(f,'C:\Users\Lara Godbersen\Documents\GitHub\Masters-thesis\Plots\Final\','beta_c', 'fontsize',17, 'figsize', [0 0 10 6.5]);

f = figure;
topoplot(data_delta_pcs.data,EEG_epoched_5.chanlocs,'colormap',viridis,'electrodes','on','emarker2', {roi_delta,'o','w',3,1})
hc=colorbar;
caxis([0.172 1.92])
xlabel(hc,'delta Power [μV^2]');
%title ('with PCS ');
set(findobj(gca,'type','patch'),'facecolor', '#F59541'); % Change [0.5, 0.5, 0.5] to your desired color
set(gca, 'FontSize', 17);
save_fig(f,'C:\Users\Lara Godbersen\Documents\GitHub\Masters-thesis\Plots\Final\','delta_pcs', 'fontsize',17, 'figsize', [0 0 10 6.5]);

f = figure;
topoplot(data_delta_c.data,EEG_epoched_5.chanlocs,'colormap',viridis,'electrodes','on','emarker2', {roi_delta,'o','w',3,1});
hc=colorbar;
caxis([0.172 1.92]);
xlabel(hc,'delta Power [μV^2]');
%title ('without PCS ');
set(findobj(gca,'type','patch'),'facecolor', '#02CAF5'); % Change [0.5, 0.5, 0.5] to your desired color
set(gca, 'FontSize', 17);
save_fig(f,'C:\Users\Lara Godbersen\Documents\GitHub\Masters-thesis\Plots\Final\','delta_c', 'fontsize',17, 'figsize', [0 0 10 6.5]);

%% 1.2 topoplots aperiodic components
data_ape_pcs = importdata('data\analysis_power\export_ape_pcs.txt');
data_ape_c = importdata('data\analysis_power\export_ape_c.txt');
data_apo_pcs = importdata('data\analysis_power\export_apo_pcs.txt');
data_apo_c = importdata('data\analysis_power\export_apo_c.txt');

% significant channel (below 5 %): 68, 86
roi_sig = [65,83];
% unter 10 % 1,2,5,6,7,35,43,65,69,71,72,78,79,80,86
roi_nearly_sig = [1,2,5,6,7,32,40,62,66,68,69,75,76,77,83];

figure;
topoplot(data_ape_pcs.data,EEG_epoched_5.chanlocs,'colormap',viridis,'electrodes','on','emarker2',{roi_sig,'o','w',3,1});
hc=colorbar;
caxis([0.648 1.16]);
xlabel(hc,'aperiodic exponent');
title ('with PCS ');
set(findobj(gca,'type','patch'),'facecolor', '#F59541'); % Change [0.5, 0.5, 0.5] to your desired color
set(gca, 'FontSize', 14);


figure;
topoplot(data_ape_c.data,EEG_epoched_5.chanlocs,'colormap',viridis,'electrodes','on','emarker2',{roi_sig,'o','w',3,1});
hc=colorbar;
caxis([0.605 1.16]);
xlabel(hc,'aperiodic exponent');
title ('without PCS ');
set(findobj(gca,'type','patch'),'facecolor', '#02CAF5'); % Change [0.5, 0.5, 0.5] to your desired color
set(gca, 'FontSize', 14);

figure;
topoplot(data_apo_pcs.data,EEG_epoched_5.chanlocs,'colormap',viridis,'electrodes','on');
hc=colorbar;
caxis([-0.630 -0.0323]);
xlabel(hc,'aperiodic offset');
title ('with PCS ');
set(findobj(gca,'type','patch'),'facecolor', '#F59541'); % Change [0.5, 0.5, 0.5] to your desired color
set(gca, 'FontSize', 16);

figure;
topoplot(data_apo_c.data,EEG_epoched_5.chanlocs,'colormap',viridis,'electrodes','on');
hc=colorbar;
caxis([-0.630 -0.0323]);
xlabel(hc,'aperiodic offset');
title ('without PCS ');
set(findobj(gca,'type','patch'),'facecolor', '#02CAF5'); % Change [0.5, 0.5, 0.5] to your desired color
set(gca, 'FontSize', 16);

   
%% 1.3 topoplot r squared
data_r_pcs = importdata('data\analysis_power\export_r_pcs.txt');
data_r_c = importdata('data\analysis_power\export_r_c.txt');


figure;
topoplot(data_r_pcs.data,EEG_epoched_5.chanlocs,'colormap',cividis,'electrodes','on');
hc=colorbar;
caxis([0.833 0.959]);
xlabel(hc,'r squared');
title ('mit PCS ');

figure;
topoplot(data_r_c.data,EEG_epoched_5.chanlocs,'colormap',cividis,'electrodes','on');
hc=colorbar;
caxis([0.833 0.959]);
xlabel(hc,'R^2');
%set(findobj(gca,'type','patch'),'facecolor', '#FFFFFF'); % Change [0.5, 0.5, 0.5] to your desired color


%% 2. plot the whole spectrum!
% for this I need the groups again
% load the age matched data
data_behav = readtable("C:\Users\Lara Godbersen\Documents\GitHub\Masters-thesis\data\PuG\matched_participants_conn.tsv", "FileType","text",'Delimiter', '\t');

% find out who is in which group
withPCS = table2array(data_behav(1:23,1));
withoutPCS = table2array(data_behav(24:46,1));

% trying to calculate a mean of the matrices
M_pcs = zeros(126,149);
M_c = zeros(126,149);
M_rest = zeros(126,149);

% load the results of the power analysis: original, oscillatory, fractal

for i = 1:length(oscillatory)
    participant_id =  oscillatory{i}.id;
    
    if ismember(participant_id,withPCS) == 1
   M_pcs =  M_pcs + oscillatory{i}.powspctrm;
    elseif ismember(participant_id,withoutPCS) == 1
   M_c = M_c + oscillatory{i}.powspctrm;
    else
        M_rest = M_rest + oscillatory{i}.powspctrm;
    end
end

avg_power_pcs = M_pcs/23;
avg_power_c = M_c/23;

v = 1;
 hold on;
    plot((original{1,1}.freq),avg_power_pcs,'k');
    plot((original{1,1}.freq),avg_power_c,'b');
    plot((fractal{v,1}.freq), mean(fractal{v,1}.powspctrm),'b');
    plot((oscillatory{v,1}.freq), mean(oscillatory{v,1}.powspctrm),'r');

    plot((original{1,1}.freq),avg_power_c,'k');
    
% now average over ROIs!
frontal_roi_c = mean(avg_power_c(roi_delta,:));
central_roi_c = mean(avg_power_c(roi_beta,:));

frontal_roi_pcs = mean(avg_power_pcs(roi_delta,:));
central_roi_pcs = mean(avg_power_pcs(roi_beta,:));

hold on;
    plot((original{1,1}.freq),frontal_roi_c,'k');% I used 'original' here just to get the freq, the data shown is from oscillatory!
    plot((original{1,1}.freq),frontal_roi_pcs,'b');
    
hold on;
    plot((original{1,1}.freq),central_roi_c,'k');
    plot((original{1,1}.freq),central_roi_pcs,'b');
    
% plot the spectral parameterization as an example for the Thesis
v = 35;
 hold on;
    plot(log((original{v,1}.freq)),log(mean(original{v,1}.powspctrm)),'k');
    plot(log((fractal{v,1}.freq)), log(mean(fractal{v,1}.powspctrm)),'b');
    plot(log((oscillatory{v,1}.freq)), log(mean(oscillatory{v,1}.powspctrm)),'r');
    
% Example RGB color code (replace with your desired color)
rgbColor = '#BA4BF5'; % This is a teal color, for example
% % Compute the logarithm of the frequency and power spectrum
logFreq = log(original{v,1}.freq);
logPower = log(mean(original{v,1}.powspctrm));

% Plot the data with the desired modifications
plot(logFreq, logPower, 'Color', rgbColor, 'LineWidth', 2);
% Add axis labels
xlabel('log(frequency)');
ylabel('log(power)');
% Set the font size for the axis ticks
set(gca, 'FontSize', 14);

% Example RGB color code (replace with your desired color)
rgbColor = '#6547F5'; % This is a teal color, for example
logFreq = log(fractal{v,1}.freq);
logPower = log(mean(fractal{v,1}.powspctrm));
% Plot the data with the desired modifications
plot(logFreq, logPower, 'Color', rgbColor, 'LineWidth', 2);
% Add axis labels
xlabel('log(frequency)');
ylabel('log(power)');
% Set the font size for the axis ticks
set(gca, 'FontSize', 14);

% Example RGB color code (replace with your desired color)
rgbColor = '#F54B46'; % This is a teal color, for example
logFreq = log(oscillatory{v,1}.freq);
logPower = log(mean(oscillatory{v,1}.powspctrm));
% Plot the data with the desired modifications
plot(logFreq, logPower, 'Color', rgbColor, 'LineWidth', 2);
% Add axis labels
xlabel('log(frequency)');
ylabel('log(power)');
% Set the font size for the axis ticks
set(gca, 'FontSize', 14);

%% 3. permutation test
% add group to every individual data set
for i = 1:length(oscillatory)
    participant_id =  oscillatory{i}.id;
    
    if ismember(participant_id,withPCS) == 1
   oscillatory{i}.group = 'withPCS';
    elseif ismember(participant_id,withoutPCS) == 1
   oscillatory{i}.group = 'withoutPCS';
    else
       oscillatory{i}.group = 'irrelevant';
    end
end

% Initialize empty arrays for P and W groups
    oscillatory_P = struct([]);
    oscillatory_W = struct([]);
% Loop through the structure array and separate based on group
for i = 1:length(oscillatory)
    if strcmp(oscillatory{i}.group, 'withPCS')
        oscillatory_P = [oscillatory_P, oscillatory{i}];
    elseif strcmp(oscillatory{i}.group, 'withoutPCS')
        oscillatory_W = [oscillatory_W, oscillatory{i}];
    end
end

% Initialize empty cell arrays for P and W groups
oscillatory_P = {};
oscillatory_W = {};

% Loop through the cell array and separate based on group
for i = 1:length(oscillatory)
    if strcmp(oscillatory{i}.group, 'withPCS')
        oscillatory_P = [oscillatory_P, oscillatory(i)];
    elseif strcmp(oscillatory{i}.group, 'withoutPCS')
        oscillatory_W = [oscillatory_W, oscillatory(i)];
    end
end

% load the layout
elec = fullfile(proj_dir,'Src\BC-128-pass-lay.mat');  
cfg = [];
cfg.elec = elec;
layout = ft_prepare_layout(cfg);
    
% prepare neighbours
channels_to_exclude = {'31', '32'};
all_channels = layout.label; % Get all channel labels from the layout
channels_to_include = setdiff(all_channels, channels_to_exclude);

    
cfg = []; 
cfg.method = 'distance'; % how should the neighbors be selected?
cfg.channel = channels_to_include;
%cfg.neighbourdist = 3; % I have no Idea what range this has, just make sure, that you get meaningful neighbors
cfg.elec = elec; % Here we need the 3d-positions!
        
neigh = ft_prepare_neighbours(cfg); % between 5 and 10 neighbors is a good value, always good to check!

% now the actual permutation test with clusters
cfg = [];
%cfg.latency          = 'all';
cfg.frequency        = [0.6 4];% [0.5 4] for delta [14 30] for beta
cfg.avgoverfreq      = 'no';
cfg.method           = 'montecarlo';
cfg.statistic        = 'ft_statfun_indepsamplesT';
cfg.parameter        = 'powspctrm';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan        = 2;
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.05;
cfg.numrandomization = 1000;
% prepare_neighbours determines what sensors may form clusters
cfg.neighbours       = neigh;

design = zeros(1,length(oscillatory_W) + length(oscillatory_P));
design(1,1:length(oscillatory_W)) = 1;
design(1,(length(oscillatory_W)+1):(length(oscillatory_W)+length(oscillatory_P))) = 2;

cfg.design           = design;
cfg.ivar             = 1;

[stat_delta] = ft_freqstatistics(cfg, oscillatory_P{:},oscillatory_W{:});
        
% get the topography of stat
plot_stat_delta = mean(stat_delta.stat,2);

f = figure;
topoplot(plot_stat_delta,EEG_epoched_5.chanlocs,'colormap',viridis,'electrodes','on','emarker2', {roi_delta,'o','w',3,1});
hc=colorbar;
caxis([-1.3711 1.3711]);
xlabel(hc,'\it{t} - value');
%title ('without PCS ');
set(gca, 'FontSize', 15);
set(findobj(gca,'type','patch'),'facecolor', '#FFFFFF'); % Change [0.5, 0.5, 0.5] to your desired color
save_fig(f,'C:\Users\Lara Godbersen\Documents\GitHub\Masters-thesis\Plots\Final\','delta_perm', 'fontsize',17, 'figsize', [0 0 10 6.5]);


% now beta
% now the actual permutation test with clusters
cfg = [];
%cfg.latency          = 'all';
cfg.frequency        = [14 30];% [0.5 4] for delta [14 30] for beta
cfg.avgoverfreq      = 'no';
cfg.method           = 'montecarlo';
cfg.statistic        = 'ft_statfun_indepsamplesT';
cfg.parameter        = 'powspctrm';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan        = 2;
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.numrandomization = 1000;
% prepare_neighbours determines what sensors may form clusters
cfg.neighbours       = neigh;

design = zeros(1,length(oscillatory_W) + length(oscillatory_P));
design(1,1:length(oscillatory_W)) = 1;
design(1,(length(oscillatory_W)+1):(length(oscillatory_W)+length(oscillatory_P))) = 2;

cfg.design           = design;
cfg.ivar             = 1;

[stat_beta] = ft_freqstatistics(cfg, oscillatory_P{:},oscillatory_W{:});

plot_stat_beta = mean(stat_beta.stat,2);
        
% get the topography of stat
f = figure;
topoplot(plot_stat_beta,EEG_epoched_5.chanlocs,'colormap',viridis,'electrodes','on','emarker2', {roi_beta,'o','w',3,1});
hc=colorbar;
caxis([-1.3711 1.3711]);
xlabel(hc,'\it{t} - value');
%title ('without PCS ');
set(gca, 'FontSize', 15);
set(findobj(gca,'type','patch'),'facecolor', '#FFFFFF'); % Change [0.5, 0.5, 0.5] to your desired color
save_fig(f,'C:\Users\Lara Godbersen\Documents\GitHub\Masters-thesis\Plots\Final\','beta_perm', 'fontsize',17, 'figsize', [0 0 10 6.5]);
