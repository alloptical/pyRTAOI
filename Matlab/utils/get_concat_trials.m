function [conc_data,trial_idx,trial_count] = get_concat_trials(data,fds)
conc_data = [];
trial_count = 0;
for f = 1:numel(fds)
    this_fd = fds{f};
    this_data = data.(this_fd)';
    conc_data = [conc_data;this_data(:)];
    this_num_trials = size(this_data,2);
    trial_idx.(this_fd) = trial_count+ [1:this_num_trials];
    trial_count = trial_count + this_num_trials;
end
end

