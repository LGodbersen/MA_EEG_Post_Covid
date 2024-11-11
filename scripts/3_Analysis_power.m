%% Power Analysis Script for my Master's thesis
% 0. Preliminaries
% 1. Spectral Parameterization
% 2. create table with all relevant values
% 3. adding group membership and other behavioral data
% 4. save combined table as .csv
% 5. fitting knee
% 6. create a big table with the resulst of the python MATLAB wrapper
% 7. add group membership & data_behav
% 8. save

%% 0.Preliminaries
clear
close all
clc

% initialize fieldtrip
ft_defaults;
% initialize EEGLAB
eeglab;close;
% add source folder for functions, toolboxes, etc.
addpath(fullfile(pwd,"Src"))
% needs MATLABs Optimization Toolbox as well! (downloaded as add on)

% set paths
proj_dir = fullfile(pwd); % automatically get path of script location, and parent dir
indir = fullfile(proj_dir,'data\prep_power_01');% path to folder with 5s epoched preprocessed data sets
outdir = fullfile(proj_dir,'data\analysis_power'); 
indat = dir(indir); % content of that folder
indat = indat(startsWith({indat.name}, 'sub-')); % only keep folders that start with 'sub-' (i.e. the subjects)

fractal = cell(length(indat),1);
original = cell(length(indat),1);
oscillatory = cell(length(indat),1); 

for s = 1:length(indat) % loop over all subjects
%% 1. Spectral Parameterization
% is originally Python Code, but was implemented in the Brainstorm Toolbox
% but not in EEGLAB -> convert back to FieldTrip
load(fullfile(indir,indat(s).name));
data_prep = eeglab2fieldtrip(EEG_epoched_5, 'raw');% the 'raw' specification seems to be important

% https://www.fieldtriptoolbox.org/example/fooof/
% compute the fractal and original spectra
    cfg               = [];
    cfg.foi        = 0.3 :0.2: 30;% changed to 0.5 for better fit
    cfg.method        = 'mtmfft';
    cfg.taper         = 'hanning';
    cfg.output        = 'fooof_aperiodic'; % there is also fooof_peaks and just fooof
    fractal{s} = ft_freqanalysis(cfg, data_prep); % this contains info an the aperiodic component already -> I do not need fooof option additionally
    cfg.output        = 'pow';
    original{s} = ft_freqanalysis(cfg, data_prep);% contains the absolute power
    
    fractal{s}.id = extractBefore(indat(s).name,'_');% attach ID

    % subtract the fractal component from the power spectrum
    cfg               = [];
    cfg.parameter     = 'powspctrm';
    cfg.operation     = 'x2-x1';
    oscillatory{s} = ft_math(cfg, fractal{s}, original{s});% contains the relative power
    
    oscillatory{s}.id = extractBefore(indat(s).name,'_');% attach ID
end

% save oscillatory, original and fractal for Visualize_power (Script)
save(fullfile(outdir,"power_results_final_01.mat"),"oscillatory", "original", "fractal");

    %% 2. create the table with all the relevant values
    % f.ex. ID, channel name (label), aperiodic component (both offset and exponent), absolute delta,
    % relative delta, absolute beta, relative beta, ...
 
% Initialize an empty table
bigTable = table();

   for s = 1:length(oscillatory)
    % subject ID
    clear patient_id
    patient_id = oscillatory{s,1}.id; 
    
     % aperiodic component (needs to be computed before the loop ->
     % separate the two values)
    clear aperiodic
    aperiodic = cell(126,1);
    for row = 1:126
        aperiodic{row,1} = fractal{s,1}.fooofparams(row).aperiodic_params;
    end
    
    aperiodic_offset = cell(126,1);
    aperiodic_exponent = cell(126,1);
    
    for row = 1:126
        aperiodic_offset{row,1} = aperiodic{row}(1);
    end
    
    for row = 1:126
        aperiodic_exponent{row,1} = aperiodic{row}(2);
    end
    
    % getting the r2 as a measure for fooof fit into the table
    for row = 1:126
        r_squared{row,1} = fractal{s,1}.fooofparams(row).r_squared;
    end
 
    % create empty cell
    cell_info = cell(126,13); 
    
for row = 1:126
   for col = 1
      cell_info{row,col} = patient_id;% VPCode
   end 
   for col = 2
       cell_info{row,col} = oscillatory{s,1}.label(row);% channel labels
   end
   for col = 3
       cell_info{row,col} = aperiodic_offset{row};% first value in fooofparams.aperiodic_parameters
   end
   for col = 4
       cell_info{row,col} = aperiodic_exponent{row};% second value in fooofparams.aperiodic_parameters
   end
   for col = 5
       cell_info{row,col} = sum(original{s,1}.powspctrm(row,4:19),2);% relevant frequencies for absoulte delta power
   end
   for col = 6
       cell_info{row,col} = sum(oscillatory{s,1}.powspctrm(row,4:19),2);% relevant frequencies for relative delta power
   end
   for col = 7
       cell_info{row,col} = sum(original{s,1}.powspctrm(row,70:149),2);% relevant frequencies for absoulte beta power
   end
   for col = 8
       cell_info{row,col} = sum(oscillatory{s,1}.powspctrm(row,70:149),2);% relevant frequencies for relative beta power
   end 
   for col = 9
       cell_info{row,col} = r_squared{row};% first value in fooofparams.aperiodic_parameters
   end
   for col = 10
       cell_info{row,col} = sum(oscillatory{s,1}.powspctrm(row,20:39),2);% relevant frequencies for relative theta power
   end
   for col = 11
       cell_info{row,col} = sum(oscillatory{s,1}.powspctrm(row,40:69),2);% relevant frequencies for relative alpha power
   end
   for col = 12
       cell_info{row,col} = sum(oscillatory{s,1}.powspctrm(row,70:99),2);% relevant frequencies for relative low beta power
   end
   for col = 13
       cell_info{row,col} = sum(oscillatory{s,1}.powspctrm(row,100:149),2);% relevant frequencies for relative high beta power
   end
end 

% create the table names
VarNames = ["participant_id" "channel" "aperiodic_offset" "aperiodic_exponent" "abs_delta" "rel_delta" "abs_beta" "rel_beta" "r_squared" "rel_theta" "rel_alpha" "rel_beta1" "rel_beta2"];

% gather data in a temporary table
currTable = table(cell_info(:,1),cell_info(:,2),cell_info(:,3),cell_info(:,4),cell_info(:,5),cell_info(:,6),cell_info(:,7),cell_info(:,8),cell_info(:,9), cell_info(:,10),cell_info(:,11),cell_info(:,12),cell_info(:,13), 'VariableNames',VarNames);
    
% Vertically concatenate the new table to the big table
bigTable = vertcat(bigTable, currTable);

   end

   %% 3. adding group membership and other behavioral data
data_behav = readtable("C:\Users\Lara Godbersen\Documents\GitHub\Masters-thesis\data\PuG\participants.tsv", "FileType","text",'Delimiter', '\t');

% combine the tables
bigTable_combined = outerjoin(bigTable, data_behav,'MergeKeys',true);

   %% 4. saving big table as .csv
% export table as csv -> be able to import it in R
csvFile = 'table_power_final_01.csv';
writetable(bigTable_combined, fullfile(outdir,csvFile));

%% 5. fitting knee 
% adding the MATLAB Python Fooof Wrapper so that I can fit 'knee' mode
% https://irenevigueguix.wordpress.com/2020/03/25/loading-python-into-matlab/
% https://github.com/fooof-tools/fooof_mat/issues/29
for v = 1:length(original)
    m{v} = mean(original{v}.powspctrm,[3],'omitnan'); %average powerspectrum over time
end

for s = 1:length(indat)
load(fullfile(indir,indat(s).name));
data_prep{s} = eeglab2fieldtrip(EEG_epoched_5, 'raw');
end

% make sure that python is loaded: type in pyversion
% make sure that the Src folder is added to the path!

clear v
for v = 1:length(original)
    settings = [];
    settings.aperiodic_mode = 'knee';
    f_range = [0.5 30]; %fitting range
    return_model = 1;
    freqs{v} = original{v}.freq;
    
    for c = 1:length(data_prep{v}.label)
        power_spectrum{v}(c,:) = m{v}(c,:);
        try
            fooof_results{v}(c,:) = fooof(freqs{v}, power_spectrum{v}(c,:) , f_range ,settings , return_model);
        catch ME
            continue
        end
    end

end

save(fullfile(outdir,"fooof_results_object_01.mat"),"fooof_results");

%% 6. put together a big table with all the fooof results
bigTable_fooof = table();

   for s = 1:length(oscillatory) % das hier über oscillatory machen
       
    % subject ID
    clear patient_id
    patient_id = oscillatory{s,1}.id; % und hier dann die ID aus fractal/oscillatorx
    
     % aperiodic component (needs to be computed before the loop ->
     % separate the two values)
    clear aperiodic
    aperiodic = cell(126,1);
    for row = 1:126
        aperiodic{row,1} = fooof_results{1,s}(row).aperiodic_params;
    end
    
    aperiodic_offset = cell(126,1);
    aperiodic_exponent = cell(126,1);
    aperiodic_knee = cell(126,1);
    
        for row = 1:126
    % Check if aperiodic is empty or not
            if isempty(aperiodic{row})
         aperiodic_offset{row} = NaN(1, 1); % Fill with NaN
            end
    
    % Check if aperiodic is empty or not
            if ~isempty(aperiodic{row})
            aperiodic_offset{row, 1} = aperiodic{row}(1);
            end
        end
        
    for row = 1:126
    % Check if aperiodict is empty or not
        if isempty(aperiodic{row})
        aperiodic_exponent{row} = NaN(1, 1); % Fill with NaN
        end
    
    % Check if aperiodic is empty or not
        if ~isempty(aperiodic{row})
        aperiodic_exponent{row, 1} = aperiodic{row}(3);
        end
    end
    
    for row = 1:126
    % Check if aperiodic is empty or not
        if isempty(aperiodic{row})
        aperiodic_knee{row} = NaN(1,1); % Fill with NaN
        end
    
    % Check if aperiodic is empty or not
        if ~isempty(aperiodic{row})
        aperiodic_knee{row, 1} = aperiodic{row}(2);
        end
    end
    
     % create powerspectrum without aperiodic component 
     for row = 1:126
    % Check if aperiodic is empty or not
        if isempty(fooof_results{1,s}(row).power_spectrum)
        spectrum_wo_ap{row,1} = NaN(1, 148); % Fill with NaN
        end
        if ~isempty(fooof_results{1,s}(row).power_spectrum)
        spectrum_wo_ap{row,1} = fooof_results{1,s}(row).power_spectrum - fooof_results{1,s}(row).ap_fit;
        end
     end
      % Power absolute
      for row = 1:126
        if isempty(fooof_results{1,s}(row).power_spectrum)
            pow_abs{row,1} = NaN(1, 148); % Fill with NaN
        end
        if ~isempty(fooof_results{1,s}(row).power_spectrum)
        pow_abs{row,1}  = fooof_results{1,s}(row).power_spectrum;
        end
      end            
                
    % create empty cell
    cell_info = cell(126,9); 
    
for row = 1:126
   for col = 1
      cell_info{row,col} = patient_id;% VPCode
   end 
   for col = 2
       cell_info{row,col} = oscillatory{s,1}.label(row);% channel labels
   end
   for col = 3
       cell_info{row,col} = aperiodic_offset{row};% first value in fooofparams.aperiodic_parameters
   end
   for col = 4
       cell_info{row,col} = aperiodic_exponent{row};% second value in fooofparams.aperiodic_parameters
   end
   for col = 5
       cell_info{row,col} = sum(pow_abs{row,1}(2:19),'omitnan');% relevant frequencies for absoulte delta power
   end
   for col = 6
       cell_info{row,col} = sum(spectrum_wo_ap{row,1}(2:19));% relevant frequencies for relative delta power
   end
   for col = 7
       cell_info{row,col} = sum(pow_abs{row,1}(70:148));% relevant frequencies for absoulte beta power
   end
   for col = 8
       cell_info{row,col} = sum(spectrum_wo_ap{row,1}(70:148));% relevant frequencies for relative beta power
   end 
   for col = 9
       cell_info{row,col} = aperiodic_knee{row};
   end 
end 

% create the table names
VarNames = ["participant_id" "channel" "aperiodic_offset" "aperiodic_exponent" "abs_delta" "rel_delta" "abs_beta" "rel_beta" "aperiodic_knee"];

% gather data in a temporary table
currTable = table(cell_info(:,1),cell_info(:,2),cell_info(:,3),cell_info(:,4),cell_info(:,5),cell_info(:,6),cell_info(:,7),cell_info(:,8),cell_info(:,9), 'VariableNames',VarNames);
    
% Vertically concatenate the new table to the big table
bigTable_fooof = vertcat(bigTable_fooof, currTable);

   end

%% 7.   adding behavioral data
%Convert your table to an array of doubles
dataAsArray = table2array(bigTable_fooof);

test =  rmmissing(bigTable_fooof);
bigTable_fooof = standardizeMissing(bigTable_fooof, {'NaN','NA', 'na', 'N/A', 'n/a', ''});
test =  rmmissing(bigTable_fooof);

% Check for NaN values in the array
nanRows = any(isnan(dataAsArray), 2);
bigTable_fooof = bigTable_fooof(~missingRows, :);

bigTable_combined_fooof = outerjoin(bigTable_fooof, data_behav,'MergeKeys',true);

   %% 8. saving big table as .csv
% export table as csv -> be able to import it in R
csvFile = 'table_power_knee_fit_01.csv';
writetable(bigTable_combined_fooof, fullfile(pwd,'data','analysis_power',csvFile));
