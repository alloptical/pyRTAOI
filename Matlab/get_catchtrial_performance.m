function [pc_lick,lick1,lick2,num_trials] = get_catchtrial_performance(trial_indices)
% fraction lick, lick left, lick right
% catch trial, photo1, photo2
pc_lick = struct();
lick1 = struct();
lick2 = struct();
num_trials = struct();

% non photo trials
this_num_trials = numel(trial_indices.stim_5_nonphotostim);
this_num_lick1 = numel([trial_indices.stim_3_var_5_correct,trial_indices.stim_4_var_5_incorrect]);
this_num_lick2 = numel([trial_indices.stim_3_var_5_incorrect,trial_indices.stim_4_var_5_correct]);
this_pc_lick = (this_num_lick1 + this_num_lick2)/this_num_trials;
this_lick1 = this_num_lick1/this_num_trials;
this_lick2 = this_num_lick2/this_num_trials;
num_trials.('catch_nonphoto') = this_num_trials;
pc_lick.catch_nonphoto = this_pc_lick;
lick1.catch_nonphoto = this_lick1;
lick2.catch_nonphoto = this_lick2;

% photo 1 trials
this_num_trials = numel(trial_indices.stim_5_photostim_1);
this_num_lick1 = numel([trial_indices.stim_3_var_3_correct,trial_indices.stim_4_var_3_incorrect]);
this_num_lick2 = numel([trial_indices.stim_3_var_3_incorrect,trial_indices.stim_4_var_3_correct]);
this_pc_lick = (this_num_lick1 + this_num_lick2)/this_num_trials;
this_lick1 = this_num_lick1/this_num_trials;
this_lick2 = this_num_lick2/this_num_trials;
num_trials.('catch_photo1') = this_num_trials;
pc_lick.catch_photo1 = this_pc_lick;
lick1.catch_photo1 = this_lick1;
lick2.catch_photo1 = this_lick2;


% photo 2 trials
this_num_trials = numel(trial_indices.stim_5_photostim_2);
this_num_lick1 = numel([trial_indices.stim_3_var_4_correct,trial_indices.stim_4_var_4_incorrect]);
this_num_lick2 = numel([trial_indices.stim_3_var_4_incorrect,trial_indices.stim_4_var_4_correct]);
this_pc_lick = (this_num_lick1 + this_num_lick2)/this_num_trials;
this_lick1 = this_num_lick1/this_num_trials;
this_lick2 = this_num_lick2/this_num_trials;
num_trials.('catch_photo2') = this_num_trials;
pc_lick.catch_photo2 = this_pc_lick;
lick1.catch_photo2 = this_lick1;
lick2.catch_photo2 = this_lick2;

end

