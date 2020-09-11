function [overall,tex1,tex2,num_trials] = get_2afc_performance(trial_indices)
overall = struct();
tex1 = struct();
tex2 = struct();
num_trials = struct();

num_correct2 = numel(trial_indices.stim_2_var_2_correct);
num_incorrect2 =  numel(trial_indices.stim_2_var_2_incorrect);
num_correct1 = numel(trial_indices.stim_1_var_1_correct);
num_incorrect1 = numel(trial_indices.stim_1_var_1_incorrect);

overall.all = (num_correct1+num_correct2)/(num_correct1+num_correct2+num_incorrect1+num_incorrect2);
tex1.all = num_correct1/(num_correct1+num_incorrect1);
tex2.all = num_correct2/(num_correct2+num_incorrect2);

try
num_correct2 = numel(intersect(trial_indices.stim_2_var_2_correct,trial_indices.stim_2_var_2_photostim));
num_incorrect2 = numel(intersect(trial_indices.stim_2_var_2_incorrect,trial_indices.stim_2_var_2_photostim));
num_correct1 = numel(intersect(trial_indices.stim_1_var_1_correct,trial_indices.stim_1_var_1_photostim));
num_incorrect1 = numel(intersect(trial_indices.stim_1_var_1_incorrect,trial_indices.stim_1_var_1_photostim));
num_trials.('stim_1_photo') = numel(trial_indices.stim_1_var_1_photostim);
num_trials.('stim_2_photo') =  numel(trial_indices.stim_2_var_2_photostim);

overall.photo = (num_correct1+num_correct2)/(num_correct1+num_correct2+num_incorrect1+num_incorrect2);
tex1.photo = num_correct1/(num_correct1+num_incorrect1);
tex2.photo = num_correct2/(num_correct2+num_incorrect2);

num_correct2 = numel(intersect(trial_indices.stim_2_var_2_correct,trial_indices.stim_2_var_2_nonphotostim));
num_incorrect2 =  numel(intersect(trial_indices.stim_2_var_2_incorrect,trial_indices.stim_2_var_2_nonphotostim));
num_correct1 = numel(intersect(trial_indices.stim_1_var_1_correct,trial_indices.stim_1_var_1_nonphotostim));
num_incorrect1 = numel(intersect(trial_indices.stim_1_var_1_incorrect,trial_indices.stim_1_var_1_nonphotostim));
num_trials.('stim_1_nonphoto') = numel(trial_indices.stim_1_var_1_nonphotostim);
num_trials.('stim_2_nonphoto') =  numel(trial_indices.stim_2_var_2_nonphotostim);

overall.nonphoto = (num_correct1+num_correct2)/(num_correct1+num_correct2+num_incorrect1+num_incorrect2);
tex1.nonphoto = num_correct1/(num_correct1+num_incorrect1);
tex2.nonphoto = num_correct2/(num_correct2+num_incorrect2);

catch
    warning('no phototrial found')
end

% dummy photo (control)
try
num_correct2 = numel(intersect(trial_indices.stim_2_var_2_correct,trial_indices.stim_2_var_2_dummyphotostim));
num_incorrect2 =  numel(intersect(trial_indices.stim_2_var_2_incorrect,trial_indices.stim_2_var_2_dummyphotostim));
num_correct1 = numel(intersect(trial_indices.stim_1_var_1_correct,trial_indices.stim_1_var_1_dummyphotostim));
num_incorrect1 = numel(intersect(trial_indices.stim_1_var_1_incorrect,trial_indices.stim_1_var_1_dummyphotostim));
num_trials.('stim_1_dummyphoto') = numel(trial_indices.stim_1_var_1_dummyphotostim);
num_trials.('stim_2_dummyphoto') =  numel(trial_indices.stim_2_var_2_dummyphotostim);


overall.dummyphoto = (num_correct1+num_correct2)/(num_correct1+num_correct2+num_incorrect1+num_incorrect2);
tex1.dummyphoto = num_correct1/(num_correct1+num_incorrect1);
tex2.dummyphoto = num_correct2/(num_correct2+num_incorrect2);catch
    warning('no dummy phototrial found')
end

% oppo photo (stimulate the other target ensemble when predicted correct)
try
num_correct2 = numel(intersect(trial_indices.stim_2_var_2_correct,trial_indices.stim_2_var_2_oppophotostim));
num_incorrect2 =  numel(intersect(trial_indices.stim_2_var_2_incorrect,trial_indices.stim_2_var_2_oppophotostim));
num_correct1 = numel(intersect(trial_indices.stim_1_var_1_correct,trial_indices.stim_1_var_1_oppophotostim));
num_incorrect1 = numel(intersect(trial_indices.stim_1_var_1_incorrect,trial_indices.stim_1_var_1_oppophotostim));
num_trials.('stim_1_oppophoto') = numel(trial_indices.stim_1_var_1_oppophotostim);
num_trials.('stim_2_oppophoto') =  numel(trial_indices.stim_2_var_2_oppophotostim);


overall.oppophoto = (num_correct1+num_correct2)/(num_correct1+num_correct2+num_incorrect1+num_incorrect2);
tex1.oppophoto = num_correct1/(num_correct1+num_incorrect1);
tex2.oppophoto = num_correct2/(num_correct2+num_incorrect2);catch
    warning('no dummy phototrial found')
end



end

