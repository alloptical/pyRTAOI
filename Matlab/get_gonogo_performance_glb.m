function [dp,hr,fa,num_trials,correct_rate] = get_gonogo_performance_glb(trial_indices)
dp = struct();
hr = struct();
fa = struct();
num_trials = struct();

num_hit = numel(trial_indices.stim_2_nonphoto_correct)+ numel(trial_indices.stim_2_photo_correct)+numel(trial_indices.stim_2_dummyphoto_correct);
num_miss = numel(trial_indices.stim_2_nonphoto_incorrect)+ numel(trial_indices.stim_2_photo_incorrect)+numel(trial_indices.stim_2_dummyphoto_incorrect);
num_cr = numel(trial_indices.stim_1_nonphoto_correct)+ numel(trial_indices.stim_1_photo_correct)+numel(trial_indices.stim_1_dummyphoto_correct);
num_fa = numel(trial_indices.stim_1_nonphoto_incorrect)+ numel(trial_indices.stim_1_photo_incorrect)+numel(trial_indices.stim_1_dummyphoto_incorrect);

[dp.all,hr.all,fa.all] = get_dprime(num_hit,num_miss,num_fa,num_cr);
correct_rate.all = (num_hit+num_cr)/(num_hit+num_cr+num_fa+num_miss);
try
num_hit = numel(trial_indices.stim_2_photo_correct);
num_miss =  numel(trial_indices.stim_2_photo_incorrect);
num_cr =  numel(trial_indices.stim_1_photo_correct);
num_fa =  numel(trial_indices.stim_1_photo_incorrect);
num_trials.('stim_1_photo') = num_cr+num_fa;
num_trials.('stim_2_photo') = num_hit+num_miss;
correct_rate.photo = (num_hit+num_cr)/(num_hit+num_cr+num_fa+num_miss);

[dp.photo,hr.photo,fa.photo] = get_dprime(num_hit,num_miss,num_fa,num_cr);

num_hit = numel(trial_indices.stim_2_nonphoto_correct);
num_miss =  numel(trial_indices.stim_2_nonphoto_incorrect);
num_cr = numel(trial_indices.stim_1_nonphoto_correct);
num_fa = numel(trial_indices.stim_1_nonphoto_incorrect);
num_trials.('stim_1_nonphoto') = num_cr+num_fa;
num_trials.('stim_2_nonphoto') =  num_hit+num_miss;
correct_rate.nonphoto = (num_hit+num_cr)/(num_hit+num_cr+num_fa+num_miss);

[dp.nonphoto,hr.nonphoto,fa.nonphoto,correct_rate.nonphoto] = get_dprime(num_hit,num_miss,num_fa,num_cr);

num_hit = numel(trial_indices.stim_2_dummyphoto_correct);
num_miss =  numel(trial_indices.stim_2_dummyphoto_incorrect);
num_cr = numel(trial_indices.stim_1_dummyphoto_correct);
num_fa = numel(trial_indices.stim_1_dummyphoto_incorrect);
num_trials.('stim_1_dummyphoto') = num_cr+num_fa;
num_trials.('stim_2_dummyphoto') =  num_hit+num_miss;
correct_rate.dummyphoto = (num_hit+num_cr)/(num_hit+num_cr+num_fa+num_miss);

[dp.dummyphoto,hr.dummyphoto,fa.dummyphoto,correct_rate.dummyphoto] = get_dprime(num_hit,num_miss,num_fa,num_cr);

catch
    warning('no phototrial found')
end

end

