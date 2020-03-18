function [dp,hr,fa] = get_gonogo_performance(trial_indices)
dp = struct();
hr = struct();
fa = struct();

num_hit = numel(trial_indices.stim_2_var_2_correct);
num_miss =  numel(trial_indices.stim_2_var_2_incorrect);
num_cr = numel(trial_indices.stim_1_var_1_correct);
num_fa = numel(trial_indices.stim_1_var_1_incorrect);

[dp.all,hr.all,fa.all] = get_dprime(num_hit,num_miss,num_fa,num_cr);

try

num_hit = numel(intersect(trial_indices.stim_2_var_2_correct,trial_indices.stim_2_var_2_photostim));
num_miss = numel(intersect(trial_indices.stim_2_var_2_incorrect,trial_indices.stim_2_var_2_photostim));
num_cr = numel(intersect(trial_indices.stim_1_var_1_correct,trial_indices.stim_1_var_1_photostim));
num_fa = numel(intersect(trial_indices.stim_1_var_1_incorrect,trial_indices.stim_1_var_1_photostim));

[dp.photo,hr.photo,fa.photo] = get_dprime(num_hit,num_miss,num_fa,num_cr);

num_hit = numel(intersect(trial_indices.stim_2_var_2_correct,trial_indices.stim_2_var_2_nonphotostim));
num_miss =  numel(intersect(trial_indices.stim_2_var_2_incorrect,trial_indices.stim_2_var_2_nonphotostim));
num_cr = numel(intersect(trial_indices.stim_1_var_1_correct,trial_indices.stim_1_var_1_nonphotostim));
num_fa = numel(intersect(trial_indices.stim_1_var_1_incorrect,trial_indices.stim_1_var_1_nonphotostim));

[dp.nonphoto,hr.nonphoto,fa.nonphoto] = get_dprime(num_hit,num_miss,num_fa,num_cr);

num_hit = numel(intersect(trial_indices.stim_2_var_2_correct,trial_indices.stim_2_var_2_dummyphotostim));
num_miss =  numel(intersect(trial_indices.stim_2_var_2_incorrect,trial_indices.stim_2_var_2_dummyphotostim));
num_cr = numel(intersect(trial_indices.stim_1_var_1_correct,trial_indices.stim_1_var_1_dummyphotostim));
num_fa = numel(intersect(trial_indices.stim_1_var_1_incorrect,trial_indices.stim_1_var_1_dummyphotostim));

[dp.dummyphoto,hr.dummyphoto,fa.dummyphoto] = get_dprime(num_hit,num_miss,num_fa,num_cr);
catch
    warning('no phototrial found')
end

end

