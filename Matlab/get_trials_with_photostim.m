function [photostim_trial_idx,num_photo_per_trial] = get_trials_with_photostim(this_trialon_frames,photo_stim_frames)
num_trials = length(this_trialon_frames);
num_photo_per_trial = zeros(1,num_trials);
for s = 1:num_trials-1
    this_start = this_trialon_frames(s);
    this_end = this_trialon_frames(s+1);
    num_photo_per_trial(s)=length(find(photo_stim_frames>this_start&photo_stim_frames<this_end));
end
photostim_trial_idx = find(num_photo_per_trial>0);
end

