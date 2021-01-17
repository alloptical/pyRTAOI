function [dp,hr,fa] = get_gonogo_baseline_performance(trial_indices)
dp = struct();
hr = struct();
fa = struct();

num_hit = numel(trial_indices.stim_2_correct);
num_miss =  numel(trial_indices.stim_2_incorrect);
num_cr = numel(trial_indices.stim_1_correct);
num_fa = numel(trial_indices.stim_1_incorrect);

[dp.all,hr.all,fa.all] = get_dprime(num_hit,num_miss,num_fa,num_cr);



end

