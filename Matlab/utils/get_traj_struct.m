function [traj_struct] = get_traj_struct(filt_cell_struct,W,trial_fds,varargin)
cell_mean = [];
for v = 1:numel(varargin)
    if strcmpi(varargin{v},'cell_mean')
        cell_mean = varargin{v+1};
    end
end
    % apply decoding and encoding matrices on raw data
    traces_in = []; %[trials, cells]
    traces_idx_struct = struct();
    for f = 1:numel(trial_fds)
        this_traces = [];
        for c = 1:size(filt_cell_struct,2)
            this_trace = filt_cell_struct(c).(trial_fds{f})(:);
            this_traces = [this_traces this_trace];
        end
        temp_size = size(traces_in,1);
        traces_in = [traces_in; this_traces];
        traces_idx_struct.(trial_fds{f}) = temp_size+1:size(traces_in,1);
        
    end
    if isempty(cell_mean)
    cell_mean =  mean(traces_in,1);
    end

    traces_in = traces_in - cell_mean; % normalise to mean
    
    F = traces_in*W;
    traj_struct = struct();
    num_dPC = size(W,2);
    trial_length = size(filt_cell_struct(1).(trial_fds{1}),1);
    
    for f = 1:numel(trial_fds)
        this_trial_idx = traces_idx_struct.(trial_fds{f});
        this_traces = F(this_trial_idx,:);
        this_num_trials = size(this_traces,1)/trial_length;
        this_traces_reshape = nan(this_num_trials,trial_length,num_dPC); %[trial,frames,factor]
        for t = 1:this_num_trials
            for m = 1:num_dPC
                this_traces_reshape(t,:,m) = this_traces([1:trial_length]+(t-1)*trial_length,m)';
            end
        end       
            traj_struct.(trial_fds{f}) = this_traces_reshape;
    end

end

