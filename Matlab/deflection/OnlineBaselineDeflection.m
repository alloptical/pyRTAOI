%% For baseline session

% stim 1 - left discrimination trials (10 conditions)
% stim 2 - right discrimination trials (10 conditions)
% stim 3 - threshold trials (+ nogo catch)
% stim 4 - sensory + photostim (10 conditions)
% stim 5 - motor only (no cue)
% stim 6 - cue only (no motor)
%% load Pybehavior data
try
[pb_file,pb_path] = uigetfile('*.mat','Select pybehavior data');
disp(['Loaded file :',fullfile(pb_path,pb_file)])
behavior_data =  load(fullfile(pb_path,pb_file)) ;
trials = make_trials_struct(behavior_data); % need to get variation type here
FLAG_PYBEHAV_LOADED = true;
catch
    FLAG_PYBEHAV_LOADED = false;
end

%% load CAIMAN data
[caiman_file,caiman_path] = uigetfile('*.mat','Select texture caiman data');
disp(['Loaded file :',fullfile(caiman_path,caiman_file)])
caiman_data = load(fullfile(caiman_path,caiman_file)); 

%% config save path
save_path = [caiman_path filesep 'analysis_files'];
if ~exist(save_path, 'dir')
    mkdir(save_path)
end

fig_save_path = [caiman_path filesep 'figures'];
if ~exist(fig_save_path, 'dir')
mkdir(fig_save_path)
end

%%
plot_stim_types = [3 1 1 2 2]; % these need to be identical to experiment_trial_seq
plot_var_types  = [2 3 4 5 6];
plot_psycho_curve(trials,plot_stim_types,plot_var_types)
%% decode choice 
% for only the difficult variations
trial_indices = struct(); % % get trialtype-outcome indices
