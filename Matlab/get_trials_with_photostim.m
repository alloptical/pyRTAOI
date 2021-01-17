function [photostim_trial_idx,num_photo_per_trial,first_photo_frames] = get_trials_with_photostim(this_trialon_frames,photo_stim_frames)
num_trials = length(this_trialon_frames);
num_photo_per_trial = zeros(1,num_trials);
first_photo_frames = nan(1,num_trials);

for s = 1:num_trials-1
    this_start = this_trialon_frames(s);
    this_end = this_trialon_frames(s+1);
    num_photo_per_trial(s)=length(find(photo_stim_frames>this_start&photo_stim_frames<this_end));
    if num_photo_per_trial(s)>0
    first_photo_frames(s) =photo_stim_frames( find(photo_stim_frames>this_start,1,'first'))-this_start;
    end
end
photostim_trial_idx = find(num_photo_per_trial>0);
end

% figure;
% stem(this_trialon_frames,ones(size(this_trialon_frames)));
% hold on;stem(photo_stim_frames,ones(size(photo_stim_frames)));
