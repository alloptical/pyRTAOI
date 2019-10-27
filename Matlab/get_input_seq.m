function [concat_seq,trial_idx,trial_count] = get_input_seq(cell_struct,cell_idx,fd_names,bin_size,varargin)
% concatinate trials for trajectory analysis input 
IF_MEDFILT = 0;
for v = 1:numel(varargin)
    if strcmpi(varargin{v},'IF_MEDFILT')
        IF_MEDFILT = varargin{v+1};
    end
end
concat_seq = [];
num_bins = floor(size(cell_struct(1).(fd_names{1}),1)/bin_size); % trial length after binning
for i = cell_idx
    this_seq = [];
    num_trials = 0;
    for f = 1:numel(fd_names)
        try
        X = cell_struct(i).(fd_names{f});
        this_num_trials = size(X,2);
        if IF_MEDFILT
            for t = 1:this_num_trials
                X(:,t) = medfilt1( X(:,t) ,3);
            end
        end
        if bin_size>1
            XX = cell2mat(arrayfun(@(x)mean(X((x-1)*bin_size+[1:bin_size],:,1)),1:num_bins,'UniformOutput', false)');
        else
            XX = X;
        end
        catch
            XX = [];
            this_num_trials = 0;
        end

        this_seq = [this_seq; XX(:)];
        num_trials = num_trials+this_num_trials;
    end
    concat_seq = [concat_seq,this_seq];
end

% get trial idx for cropping concatinated traces
trial_idx = struct();
trial_count = 0;
for f = 1:numel(fd_names)
    this_fd = fd_names{f}; 
    try
    this_num_trials =  size(cell_struct(i).(this_fd),2);
    catch
        this_num_trials = 0;
    end
    trial_idx.(this_fd) = trial_count+ [1:this_num_trials];
    trial_count = trial_count + this_num_trials;
end
end

