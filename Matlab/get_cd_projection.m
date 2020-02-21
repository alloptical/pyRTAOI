function [this_proj_struct] = get_cd_projection(this_cell_struct,this_cd,fds)
% fds = fields(this_cell_struct);
for f = 1:numel(fds)
    this_fd =fds{f};
    this_traces = [];
    for c = 1:size(this_cell_struct,2)
        this_trace = this_cell_struct(c).(this_fd)(:);
        this_traces = [this_traces this_trace];
    end
    
    this_proj = this_traces*this_cd';
    trial_length = size( this_cell_struct(c).(this_fd),1);
    this_num_trials =round( size(this_proj,1)/trial_length);
    this_traces_reshape = nan(this_num_trials,trial_length); %[trial,frames,factor]; in this case factor = cell
    for t = 1:this_num_trials
        this_traces_reshape(t,:) = this_proj([1:trial_length]+(t-1)*trial_length)';
    end
    
    this_proj_struct.(this_fd) = this_traces_reshape;
    
end
end

