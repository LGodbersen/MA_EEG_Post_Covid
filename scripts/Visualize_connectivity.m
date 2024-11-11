%% Script for visualising the coherence and connectivity results
% overview
% 1. mean connectivity matrix per group
    % 1.1 delta
    % 1.2 beta
% 2. plots for the overview chart 
% 3. topoplot
% 4. have a look at the modularity (exploratory)

%% 1. mean connectivity matrix per group
% 1. calculate and save all the individual connectivity matrices
% I want to plot the mean coherence matrix per group!
% what I need: loop that first creates the matrix and puts it into a
% structure for every person
% then I need the IDs in order to know who is in which group -> add a group
% identifier
% then I can average the correlation matrices
% and then plot the two
% https://de.mathworks.com/matlabcentral/answers/34395-by-element-average-of-multiple-matrices

clear
close all
clc

%initialize
ft_defaults;
eeglab;close;

% set paths
proj_dir = fullfile(pwd); % automatically get path of script location, and parent dir
indir = fullfile(proj_dir,'data\prep_connectivity_01');% path to folder with BIDS datasets
outdir = fullfile(proj_dir,'data\analysis_connectivity_icoh\01'); %  alternative: 'data\analysis_connectivity_icoh'
indat = dir(indir); % content of that folder
indat = indat(startsWith({indat.name}, 'sub-')); % only keep folders that start with 'sub-' (i.e. the subjects)

roi_delta = [21,102,11,37,72,36,46,79,45,19,109,24,91,90,80,89,92,93,20,47,10,56,25]; % frontal ROI (channel indices are not always the number of the actual electrode!)
roi_beta = [82,62,87,63,1,65,3,64,2,67,71,73,78,31,34,39,83,40,84,41,85,42,86,43,74,5,75,6,7,76,8,77,68,32,69,33,70]; % put channel indices here


for z = 1:length(indat)
    
    load(fullfile(indir,indat(z).name));
data_prep = eeglab2fieldtrip(EEG_epoched_4, 'raw');

    cfg       = [];
cfg.foi       = 0.3 :0.2: 30;% our high pass filter is at 0.5
cfg.method    = 'mtmfft';
cfg.taper     = 'dpss';
cfg.pad       = 'nextpow2';
cfg.output    = 'fourier';
cfg.tapsmofrq = 2;
freq{z}          = ft_freqanalysis(cfg, data_prep);

cfg           = [];
cfg.method    = 'coh';% coherence
cfg.complex   = 'imag';% compute imaginary part of coherence
coh{z}           = ft_connectivityanalysis(cfg, freq{z});

coh{z}.id =  extractBefore(indat(z).name,'_');

% connectivity matrix for delta
conn_matrix = coh{z}.cohspctrm;
conn_matrix_delta = mean(conn_matrix(:,:,3:16),3); 
conn_matrix_delta = abs(conn_matrix_delta);

% connectivity matrix for beta
conn_matrix_beta = mean(conn_matrix(:,:,58:122),3);
conn_matrix_beta = abs(conn_matrix_beta);

outdir = fullfile(proj_dir,'data\connectivity_matrices_delta');
    save(fullfile(outdir,[coh{z}.id + "_con_matrix_d.mat"]),"conn_matrix_delta");
    
outdir = fullfile(proj_dir,'data\connectivity_matrices_beta');
    save(fullfile(outdir,[coh{z}.id + "_con_matrix_b.mat"]),"conn_matrix_beta");

end

%% 1.1 delta
proj_dir = fullfile(pwd); % automatically get path of script location, and parent dir
indir = fullfile(proj_dir,'data\connectivity_matrices_delta');% path to folder with BIDS datasets
outdir = fullfile(proj_dir,'data\analysis_connectivity'); % alternative 'data\analysis_connectivity_icoh' if you want imaginary part of coherence
indat = dir(indir); % content of that folder
indat = indat(startsWith({indat.name}, 'sub-')); % only keep folders that start with 'sub-' (i.e. the subjects)

% load the age matched data
data_behav = readtable("C:\Users\Lara Godbersen\Documents\GitHub\Masters-thesis\data\PuG\matched_participants_conn.tsv", "FileType","text",'Delimiter', '\t');

% find out who is in which group
withPCS = table2array(data_behav(1:23,1));
withoutPCS = table2array(data_behav(24:46,1));

% trying to calculate a mean of the matrices
M_pcs = zeros(126,126);
M_c = zeros(126,126);
M_rest = zeros(126,126);

for i = 1:length(indat)
    load(fullfile(indir,indat(i).name));
    participant_id =  extractBefore(indat(i).name,'_');
    
    if ismember(participant_id,withPCS) == 1
   M_pcs =  M_pcs + abs(conn_matrix_delta);
    elseif ismember(participant_id,withoutPCS) == 1
   M_c = M_c + abs(conn_matrix_delta);
    else
        M_rest = M_rest + abs(conn_matrix_delta);
    end
end

% afterwards divide by number of people in the group
% e.g. 

avg_con_delta_pcs = M_pcs/23;
avg_con_delta_c = M_c/23;

% plot this
figure;
imagesc(avg_con_delta_pcs);
axis square; colorbar;
title('imaginary part of coherence delta, avg all subj')

% now we want to sort the channels so that neighbouring channels are next t
%o each other
%right_placement = {21, Gnd, 26, 96, 95, 59, 50, 93, 94, 123, 22, 105, 11, 40, 83, 10, 39, 49, 92, 48, 19, 112, 25, 128, 30, 122, 58, 104, 91, 82, 38, 69, 33, 75, 84, 97, 51, 113, 27, 114, 106, 12, 41, 4, 70, 66,1,65,68, 74, 9,47,18,111,121,64,120,57,103,90,81,37,3,67,2,34,76,85,98,52,115,60,124,23,107,13,42,77,71,35,72,36,73,80,46,17,110,24,127,29,119, 56, 102, 89, 8, 79, 7,6,78,5, 86,99,53,116,28,61,108,14,43,87, 44, 88, 45, 16, 109, 63, 126, 118, 55, 101, 15, 100, 54, 117, 125, 62}; 
% old: right_placement = [20, 126, 25, 93, 92, 56, 47, 90, 91, 120, 21, 102, 11, 37, 80, 10, 36, 46, 89, 45, 19, 109, 24, 125, 29, 119, 55, 101, 88, 79, 35, 66, 30, 72, 81, 94, 48, 110, 26, 111, 103, 12, 38, 4, 67, 63,1,62,65, 71, 9,44,18,108,118,61,117,54,100,87,78,34,3,64,2,31,73,82,95,49,112,57,121,22,104,13,39,74,68,32,69,33,70,77,43,17,107,23,124,28,116, 53, 99, 86, 8, 76, 7,6,75,5, 83,96,50,113,27,58,105,14,40,84, 41, 85, 42, 16, 106, 60, 123, 115, 52, 98, 15, 97, 51, 114, 122, 59]; 

right_placement = [123, 22, 105, 94, 21, 129, 26, 96, 128, 25, 112, 19, 95, 59, 50, 93, 11, 40, 83, 10, 92, 48, 49, 39,27, 113, 51, 97, 12, 98, 52, 115, 60, 114, 106, 30, 122, 58, 104, 18, 103, 57, 120, 64, 121, 111, 84, 75, 33, 69, 38, 82, 91, 9, 74, 68, 1, 66, 70, 4, 41, 76, 34, 2, 65, 3, 37, 81, 47, 90, 80, 73, 67, 71, 77, 85, 42, 5, 35, 72, 36, 8, 46, 13, 86, 78, 6, 7, 79, 89, 17, 110, 24, 127, 102, 45, 88, 44, 87, 43, 99, 107, 23, 124, 14, 100, 15, 101, 16, 53, 108, 54, 55, 109, 56, 119, 63, 118, 117, 61, 116, 28, 125, 62, 126, 29]; 

    % Initialize array
    transformed_numbers = zeros(size(right_placement));
    
    % transformation so that the numbers match with the index because f.ex.
    % the channel with number 20 is missing and also 31 and 32
    for i = 1:length(right_placement)
        if right_placement(i) <= 19
            transformed_numbers(i) = right_placement(i);
        elseif right_placement(i) > 19 && right_placement(i) <= 30
            transformed_numbers(i) = right_placement(i) - 1;
        else
            transformed_numbers(i) = right_placement(i) - 3;
        end
        % Debug output
    fprintf('Index %d: right_placement = %d, transformed_numbers = %d\n', ...
            i, right_placement(i), transformed_numbers(i));
    end


avg_con_delta_pcs_reordered = avg_con_delta_pcs(transformed_numbers, transformed_numbers);
avg_con_delta_c_reordered = avg_con_delta_c(transformed_numbers, transformed_numbers);

% Define the coordinates for the box
box_start = 0.5; % Start at 0.5 to include the first channel
box_end = 24.5;  % End at 24.5 to include the 24th channel
box_size = box_end - box_start;

box_start_tl = 24.5; % Start at 0.5 to include the first channel
box_end_tl = 35.5;  % End at 24.5 to include the 24th channel
box_size_tl = box_end_tl - box_start_tl;

box_start_tr = 35.5; % Start at 0.5 to include the first channel
box_end_tr = 46.5;  % End at 24.5 to include the 24th channel
box_size_tr = box_end_tr - box_start_tr;

box_start_c = 46.5; % Start at 0.5 to include the first channel
box_end_c = 83.5;  % End at 24.5 to include the 24th channel
box_size_c = box_end_c - box_start_c;

box_start_p = 83.5; % Start at 0.5 to include the first channel
box_end_p = 109.5;  % End at 24.5 to include the 24th channel
box_size_p = box_end_p - box_start_p;

box_start_o = 109.5; % Start at 0.5 to include the first channel
box_end_o = 126;  % End at 24.5 to include the 24th channel
box_size_o = box_end_o - box_start_o;


% start the figure
figure;
imagesc(avg_con_delta_pcs_reordered);
caxis([0 0.046]);
hc = colorbar;
axis square;
colormap(viridis);
xlabel(hc,'imag(coh) delta');
ylabel('channels');
title('with PCS')
set(gca, 'FontSize', 16);

% Draw the rectangle
hold on;
rectangle('Position', [box_start, box_start, box_size, box_size], 'EdgeColor', 'r', 'LineWidth', 1.5);
hold off;

hold on;
rectangle('Position', [box_start_tl, box_start_tl, box_size_tl, box_size_tl], 'EdgeColor', 'w', 'LineWidth', 1.5);
hold off;

hold on;
rectangle('Position', [box_start_tr, box_start_tr, box_size_tr, box_size_tr], 'EdgeColor', 'w', 'LineWidth', 1.5);
hold off;

hold on;
rectangle('Position', [box_start_c, box_start_c, box_size_c, box_size_c], 'EdgeColor', 'w', 'LineWidth', 1.5);
hold off;

hold on;
rectangle('Position', [box_start_p, box_start_p, box_size_p, box_size_p], 'EdgeColor', 'w', 'LineWidth', 1.5);
hold off;

hold on;
rectangle('Position', [box_start_o, box_start_o, box_size_o, box_size_o], 'EdgeColor', 'w', 'LineWidth', 1.5);
hold off;

figure;
imagesc(avg_con_delta_c_reordered);
caxis([0 0.046]);
hc = colorbar;
axis square;
colormap(viridis);
xlabel(hc,'imag(coh) delta');
ylabel('channels');
title('without PCS')
set(gca, 'FontSize', 16);

% Draw the rectangle
hold on;
rectangle('Position', [box_start, box_start, box_size, box_size], 'EdgeColor', 'r', 'LineWidth', 1.5);
hold off;

hold on;
rectangle('Position', [box_start_tl, box_start_tl, box_size_tl, box_size_tl], 'EdgeColor', 'w', 'LineWidth', 1.5);
hold off;

hold on;
rectangle('Position', [box_start_tr, box_start_tr, box_size_tr, box_size_tr], 'EdgeColor', 'w', 'LineWidth', 1.5);
hold off;

hold on;
rectangle('Position', [box_start_c, box_start_c, box_size_c, box_size_c], 'EdgeColor', 'w', 'LineWidth', 1.5);
hold off;

hold on;
rectangle('Position', [box_start_p, box_start_p, box_size_p, box_size_p], 'EdgeColor', 'w', 'LineWidth', 1.5);
hold off;

hold on;
rectangle('Position', [box_start_o, box_start_o, box_size_o, box_size_o], 'EdgeColor', 'w', 'LineWidth', 1.5);
hold off;

%% 1.2  beta
indir = fullfile(proj_dir,'data\connectivity_matrices_beta');% path to folder with BIDS datasets
indat = dir(indir); % content of that folder
indat = indat(startsWith({indat.name}, 'sub-')); % only keep folders that start with 'sub-' (i.e. the subjects)

M_beta_pcs = zeros(126,126);
M_beta_c = zeros(126,126);
M_beta_rest = zeros(126,126);

for i = 1:length(indat)
    load(fullfile(indir,indat(i).name));
    participant_id =  extractBefore(indat(i).name,'_');
    
    if ismember(participant_id,withPCS) == 1
   M_beta_pcs =  M_beta_pcs + abs(conn_matrix_beta);
    elseif ismember(participant_id,withoutPCS) == 1
   M_beta_c = M_beta_c + abs(conn_matrix_beta);
    else
        M_beta_rest = M_beta_rest + abs(conn_matrix_beta);
    end
end

avg_con_beta_pcs = M_beta_pcs/23;
avg_con_beta_c = M_beta_c/23;

% take the transformed numbers sequence from before
avg_con_beta_pcs_reordered = avg_con_beta_pcs(transformed_numbers, transformed_numbers);
avg_con_beta_c_reordered = avg_con_beta_c(transformed_numbers, transformed_numbers);

figure;
imagesc(avg_con_beta_pcs_reordered);
hc = colorbar;
axis square;
caxis([0 0.046]);%0.11
colormap(viridis);
xlabel(hc,'imag(coh) beta');
ylabel('channels');
title('with PCS')
set(gca, 'FontSize', 16);

% Draw the rectangle
hold on;
rectangle('Position', [box_start, box_start, box_size, box_size], 'EdgeColor', 'w', 'LineWidth', 1.5);
hold off;

hold on;
rectangle('Position', [box_start_tl, box_start_tl, box_size_tl, box_size_tl], 'EdgeColor', 'w', 'LineWidth', 1.5);
hold off;

hold on;
rectangle('Position', [box_start_tr, box_start_tr, box_size_tr, box_size_tr], 'EdgeColor', 'w', 'LineWidth', 1.5);
hold off;

hold on;
rectangle('Position', [box_start_c, box_start_c, box_size_c, box_size_c], 'EdgeColor', 'r', 'LineWidth', 1.5);
hold off;

hold on;
rectangle('Position', [box_start_p, box_start_p, box_size_p, box_size_p], 'EdgeColor', 'w', 'LineWidth', 1.5);
hold off;

hold on;
rectangle('Position', [box_start_o, box_start_o, box_size_o, box_size_o], 'EdgeColor', 'w', 'LineWidth', 1.5);
hold off;

figure;
imagesc(avg_con_beta_c_reordered);
hc = colorbar;
axis square;
caxis([0 0.046]);
colormap(viridis);
xlabel(hc,'imag(coh) beta');
ylabel('channels');
title('without PCS')
set(gca, 'FontSize', 16);

% Draw the rectangle
hold on;
rectangle('Position', [box_start, box_start, box_size, box_size], 'EdgeColor', 'w', 'LineWidth', 1.5);
hold off;

hold on;
rectangle('Position', [box_start_tl, box_start_tl, box_size_tl, box_size_tl], 'EdgeColor', 'w', 'LineWidth', 1.5);
hold off;

hold on;
rectangle('Position', [box_start_tr, box_start_tr, box_size_tr, box_size_tr], 'EdgeColor', 'w', 'LineWidth', 1.5);
hold off;

hold on;
rectangle('Position', [box_start_c, box_start_c, box_size_c, box_size_c], 'EdgeColor', 'r', 'LineWidth', 1.5);
hold off;

hold on;
rectangle('Position', [box_start_p, box_start_p, box_size_p, box_size_p], 'EdgeColor', 'w', 'LineWidth', 1.5);
hold off;

hold on;
rectangle('Position', [box_start_o, box_start_o, box_size_o, box_size_o], 'EdgeColor', 'w', 'LineWidth', 1.5);
hold off;
      
%% 2. plots for the overview chart 
% now I can also check per threshold, if the network is still connected and
% how the visualisation would look

% http://complexity.es/school/neuroscience/
% use the connectivity matrices, threshold them differently and check

conn_matrix_beta = conn_matrix_beta(transformed_numbers, transformed_numbers);

f = figure;
imagesc(conn_matrix_beta);
hc = colorbar;
axis square;
caxis([0 0.075]);
colormap(viridis);
xlabel(hc,'imag(coh) beta of one participant');
ylabel('Channels sorted by regions');
set(gca, 'FontSize', 16);

save_fig(f,'C:\Users\Lara Godbersen\Documents\GitHub\Masters-thesis\Plots\Visualize connectivity\','original_matrix', 'fontsize',16, 'figsize', [0 0 11.5 8.5]);

histogram(conn_matrix_beta);

    sortedValues_beta_01 = sort(abs(conn_matrix_beta(:)),'descend');
    sortedValues_beta_01 = sortedValues_beta_01(~isnan(sortedValues_beta_01));
    threshold_beta_01 = sortedValues_beta_01(floor(0.1*length(sortedValues_beta_01)));
    
    sortedValues_beta_03 = sort(abs(conn_matrix_beta(:)),'descend');
    sortedValues_beta_03 = sortedValues_beta_03(~isnan(sortedValues_beta_03));
    threshold_beta_03 = sortedValues_beta_03(floor(0.3*length(sortedValues_beta_03)));
    
    sortedValues_beta_07 = sort(abs(conn_matrix_beta(:)),'descend');
    sortedValues_beta_07 = sortedValues_beta_07(~isnan(sortedValues_beta_07));
    threshold_beta_07 = sortedValues_beta_07(floor(0.7*length(sortedValues_beta_07)));

   % plot histogram to understand the distribution of the values
    f = figure;
    histogram(conn_matrix_beta);
    ylabel('Amount [n]')
    xlabel('imag(coh) beta of one participant')
    xline(threshold_beta_01,'r-',{'Threshold 10 %'})
    xline(threshold_beta_07,'r-',{'Threshold 70 %'})
    xline(threshold_beta_03,'r-',{'Threshold 30 %'})
    set(gca, 'FontSize', 14);
    %xline(threshold_beta_09,'-',{'Threshold',num2str(round(threshold_beta_09,3))})
    save_fig(f,'C:\Users\Lara Godbersen\Documents\GitHub\Masters-thesis\Plots\Visualize connectivity\','histogram_beta', 'fontsize',16, 'figsize', [0 0 11.5 8.5]);

    threshold_matrix_beta_01 = threshold_proportional(conn_matrix_beta,0.1);

    isconnected(threshold_matrix_beta_01);
    
    f = figure;
imagesc(threshold_matrix_beta_01);
hc = colorbar;
axis square;
caxis([0 0.075]);
colormap(viridis);
xlabel(hc,'imag(coh) beta of one participant');
ylabel('Channels sorted by regions');
set(gca, 'FontSize', 14);

save_fig(f,'C:\Users\Lara Godbersen\Documents\GitHub\Masters-thesis\Plots\Visualize connectivity\','10_conn_matrix', 'fontsize',16, 'figsize', [0 0 11.5 8.5]);
 
    f = figure;
    plot(graph(threshold_matrix_beta_01))
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
    save_fig(f,'C:\Users\Lara Godbersen\Documents\GitHub\Masters-thesis\Plots\Visualize connectivity\','10_graph', 'fontsize',14, 'figsize', [0 0 11.5 8.5]);

    
    % 0.7 beta
    threshold_matrix_beta_07 = threshold_proportional(conn_matrix_beta,0.7);
    
 f =     figure;
imagesc(threshold_matrix_beta_07);
hc = colorbar;
axis square;
caxis([0 0.075]);
colormap(viridis);
xlabel(hc,'imag(coh) beta of one participant');
ylabel('Channels sorted by regions');
set(gca, 'FontSize', 14);

save_fig(f,'C:\Users\Lara Godbersen\Documents\GitHub\Masters-thesis\Plots\Visualize connectivity\','70_conn_matrix', 'fontsize',16, 'figsize', [0 0 11.5 8.5]);
 
f = figure;
plot(graph(threshold_matrix_beta_07))
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
    save_fig(f,'C:\Users\Lara Godbersen\Documents\GitHub\Masters-thesis\Plots\Visualize connectivity\','70_graph', 'fontsize',14, 'figsize', [0 0 11.5 8.5]);

    % 0.3 beta
    threshold_matrix_beta_03 = threshold_proportional(conn_matrix_beta,0.3);

    
    isconnected(threshold_matrix_beta_03)
    
    f = figure;
    imagesc(threshold_matrix_beta_03);
hc = colorbar;
axis square;
caxis([0 0.075]);
colormap(viridis);
xlabel(hc,'imag(coh) beta of one participant');
ylabel('Channels sorted by regions');
set(gca, 'FontSize', 14);
save_fig(f,'C:\Users\Lara Godbersen\Documents\GitHub\Masters-thesis\Plots\Visualize connectivity\','30_conn_matrix', 'fontsize',16, 'figsize', [0 0 11.5 8.5]);

f = figure;
plot(graph(threshold_matrix_beta_03))
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
    save_fig(f,'C:\Users\Lara Godbersen\Documents\GitHub\Masters-thesis\Plots\Visualize connectivity\','30_graph', 'fontsize',14, 'figsize', [0 0 11.5 8.5]);


  adjacency_matrix_beta_03 = weight_conversion(threshold_matrix_beta_03,'binarize');
    f =   figure;
imagesc(adjacency_matrix_beta_03);
hc = colorbar;
set(hc,'Ticks',[0 1]);
axis square;
caxis([0 1]);
colormap(viridis);
xlabel(hc,'binarized imag(coh) beta');
ylabel('Channels sorted by regions');
set(gca, 'FontSize', 14);
hold on;
rectangle('Position', [box_start_c, box_start_c, box_size_c, box_size_c], 'EdgeColor', 'r', 'LineWidth', 1.5);
hold off;
save_fig(f,'C:\Users\Lara Godbersen\Documents\GitHub\Masters-thesis\Plots\Visualize connectivity\','30_adj_matrix', 'fontsize',16, 'figsize', [0 0 11.5 8.5]);

    
    plot(graph(threshold_matrix_beta_03))
    axis off
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
    
    distance_beta_03 = distance_bin(adjacency_matrix_beta_03);
    distance_beta_roi_03 = distance_beta_03(roi_beta,roi_beta,:);
    
    f =   figure;
    imagesc(distance_beta_03);
hc = colorbar;
axis square;
caxis([0 2]);
set(hc,'Ticks',[0 1 2]);
colormap(viridis);
xlabel(hc,'Distance bins');
ylabel('Channels sorted by regions');
set(gca, 'FontSize', 14);
hold on;
rectangle('Position', [box_start_c, box_start_c, box_size_c, box_size_c], 'EdgeColor', 'r', 'LineWidth', 1.5);
hold off;
save_fig(f,'C:\Users\Lara Godbersen\Documents\GitHub\Masters-thesis\Plots\Visualize connectivity\','30_distance_matrix', 'fontsize',16, 'figsize', [0 0 11.5 8.5]);

    
    % 0.3 delta
    threshold_matrix_delta_03 = threshold_proportional(conn_matrix_delta,0.3);

    isconnected(threshold_matrix_delta_03)
    imagesc(threshold_matrix_delta_03); colorbar
    
    plot(graph(threshold_matrix_delta_03))
    axis off
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
    
    % 0.4 delta
    threshold_matrix_delta_04 = threshold_proportional(conn_matrix_delta,0.4);

    isconnected(threshold_matrix_delta_04)
    imagesc(threshold_matrix_delta_04); colorbar
    
    plot(graph(threshold_matrix_delta_04, 'upper'))
    axis off
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
    
    % 0.5 delta
    threshold_matrix_delta_05 = threshold_proportional(avg_con_beta_c,0.5);

    isconnected(threshold_matrix_delta_05)
        figure;
imagesc(threshold_matrix_delta_05);
hc = colorbar;
axis square;
caxis([0 0.075]);
title('imaginary part of coherence (beta) without PCS')
colormap(viridis);
xlabel(hc,'imaginärer Teil der beta band Kohärenz');
ylabel('Elektroden (unsortiert)');
title('ohne PCS')
    
    
    plot(graph(threshold_matrix_delta_05, 'upper'))
    axis off
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
    
    adjacency_matrix_delta_05 = weight_conversion(threshold_matrix_delta_05,'binarize');
      figure;
imagesc(adjacency_matrix_delta_05);
hc = colorbar;
axis square;
caxis([0 0.075]);
title('imaginary part of coherence (beta) without PCS')
colormap(viridis);
xlabel(hc,'imaginärer Teil der beta band Kohärenz');
ylabel('Elektroden (unsortiert)');
title('ohne PCS')

    % 0.6 delta
    threshold_matrix_delta_06 = threshold_proportional(conn_matrix_delta,0.6);

    isconnected(threshold_matrix_delta_06)
    imagesc(threshold_matrix_delta_06); colorbar
    
    plot(graph(threshold_matrix_delta_06, 'upper'))
    axis off
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
    
    % 0.7 delta
    threshold_matrix_delta_07 = threshold_proportional(conn_matrix_delta,0.7);

    isconnected(threshold_matrix_delta_07)
    imagesc(threshold_matrix_delta_07); colorbar
    
    plot(graph(threshold_matrix_delta_07, 'upper'))
    axis off
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
    
    % 0.8 delta
    threshold_matrix_delta_08 = threshold_proportional(conn_matrix_delta,0.8);

    isconnected(threshold_matrix_delta_08)
    imagesc(threshold_matrix_delta_08); colorbar
    
    plot(graph(threshold_matrix_delta_08, 'upper'))
    axis off
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
    
    
    % 0.9  delta
    threshold_matrix_delta_09 = threshold_proportional(conn_matrix_delta,0.9);

    isconnected(threshold_matrix_delta_09)
    imagesc(threshold_matrix_delta_09); colorbar
    
    plot(graph(threshold_matrix_delta_09, 'lower'))
    axis off
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
    
    %
    adjacency_matrix_delta_02 = weight_conversion(threshold_matrix_delta_01,'binarize');
    
    plot(graph(adjacency_matrix_delta_02))
    axis off
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
    
    adjacency_matrix_delta_roi_01 = adjacency_matrix_delta_02(roi_delta,roi_delta,:);
    
    isconnected(adjacency_matrix_delta_roi_01)
    
    imagesc(adjacency_matrix_delta_roi_01); colorbar
    
    plot(graph(adjacency_matrix_delta_roi_01,'lower'))
    axis off
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
    
%% 3. topoplot
% 3.1 delta
% with PCS
half_con_delta = avg_con_delta_pcs/2;
topo_con_delta = mean(half_con_delta,2);

figure;
topoplot(topo_con_delta,EEG_epoched_4.chanlocs,'colormap',viridis,'electrodes','on','emarker2', {roi_delta,'o','w',3,1});
hc=colorbar;
caxis([0.0194 0.0257]);
xlabel(hc,'imag(coh) delta');
title ('with PCS ');
set(gca, 'FontSize', 15);
set(findobj(gca,'type','patch'),'facecolor', '#FFFFFF'); % Change [0.5, 0.5, 0.5] to your desired color

% without PCS
half_con_delta_c = avg_con_delta_c/2;
topo_con_delta_c = mean(half_con_delta_c,2);

figure;
topoplot(topo_con_delta_c,EEG_epoched_4.chanlocs,'colormap',viridis,'electrodes','on','emarker2', {roi_delta,'o','w',3,1});
hc=colorbar;
caxis([0.0194 0.0257]);
xlabel(hc,'imag(coh) delta');
title ('without PCS ');
set(gca, 'FontSize', 15);
set(findobj(gca,'type','patch'),'facecolor', '#FFFFFF'); % Change [0.5, 0.5, 0.5] to your desired color

% 3.2 beta
% with PCS
half_con_beta = avg_con_beta_pcs/2;
topo_con_beta = mean(half_con_beta,2);

figure;
topoplot(topo_con_beta,EEG_epoched_4.chanlocs,'colormap',viridis,'electrodes','on','emarker2', {roi_beta,'o','w',3,1});
hc=colorbar;
caxis([0.0096 0.0141]);
xlabel(hc,'imag(coh) beta');
title ('with PCS ');
set(gca, 'FontSize', 15);
set(findobj(gca,'type','patch'),'facecolor', '#FFFFFF'); % Change [0.5, 0.5, 0.5] to your desired color

% without PCS
half_con_beta_c = avg_con_beta_c/2;
topo_con_beta_c = mean(half_con_beta_c,2);

figure;
topoplot(topo_con_beta_c,EEG_epoched_4.chanlocs,'colormap',viridis,'electrodes','on','emarker2', {roi_beta,'o','w',3,1});
hc=colorbar;
caxis([0.0096 0.0141]);
xlabel(hc,'imag(coh) beta');
title ('without PCS ');
set(gca, 'FontSize', 15);
set(findobj(gca,'type','patch'),'facecolor', '#FFFFFF'); % Change [0.5, 0.5, 0.5] to your desired color

    
%% 4. have a look at the modularity (exploratory)
community_louvain(conn_matrix_delta)

% Iterative community finetuning.
        % W is the input connection matrix.
        n  = size(avg_con_delta_pcs,1);             % number of nodes
        Mod_pcs  = 1:n;                   % initial community affiliations
        Q0 = -1; Q1 = 0;            % initialize modularity values
        while Q1-Q0>1e-5;           % while modularity increases
            Q0 = Q1;                % perform community detection
            [Mod_pcs, Q1] = community_louvain(avg_con_delta_pcs, [], Mod_pcs);
        end
        
% Iterative community finetuning.
        % W is the input connection matrix.
        n  = size(avg_con_delta_c,1);             % number of nodes
        Mod_c  = 1:n;                   % initial community affiliations
        Q0 = -1; Q1 = 0;            % initialize modularity values
        while Q1-Q0>1e-5;           % while modularity increases
            Q0 = Q1;                % perform community detection
            [Mod_c, Q1] = community_louvain(avg_con_delta_c, [], Mod_c);
        end

      [h,p,ci,stats] =  ttest(Mod_pcs,Mod_c);