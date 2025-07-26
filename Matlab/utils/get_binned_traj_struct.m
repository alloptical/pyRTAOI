function [binned_struct] = get_binned_traj_struct(traj_struct,bin_size)
fd_names = fields(traj_struct);
for i = 1:numel(fd_names)
    this_fd = fd_names{i};
    this_num_frames = size(traj_struct(1).(this_fd),1);
    
    X = traj_struct.(this_fd);
    trial_dim = 1;
    frame_dim = 2;
    state_dim = 3;
    this_num_trials = size(X,trial_dim);
    this_num_frames = size(traj_struct(1).(this_fd),frame_dim);
    num_bins = round(this_num_frames/bin_size);
    % binning
    if bin_size>1
        XX = cell2mat(arrayfun(@(x)mean(X(:,(x-1)*bin_size+[1:bin_size],:),2),1:num_bins,'UniformOutput', false));
    end
    
    
    binned_struct.(this_fd)= XX;
    
end
end

