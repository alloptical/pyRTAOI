function [concat_seq,trial_idx,trial_count] = get_input_seq(cell_struct,cell_idx,fd_names,bin_size,varargin)
% concatinate trials or trial averages for trajectory analysis input 
IF_MEDFILT = 0;
trial_dim = 2;
time_dim = 1;
fd_trial_idx = [];
avg_frames = [];
for v = 1:numel(varargin)
    if strcmpi(varargin{v},'IF_MEDFILT')
        IF_MEDFILT = varargin{v+1};
    end
    if strcmpi(varargin{v},'trial_dim')
        trial_dim = varargin{v+1};
    end
    if strcmpi(varargin{v},'time_dim')
        time_dim = varargin{v+1};
    end
    if strcmpi(varargin{v},'fd_trial_idx')
        fd_trial_idx = varargin{v+1};
    end
    if strcmpi(varargin{v},'avg_frames')
        avg_frames = varargin{v+1};
    end
end
concat_seq = [];
num_bins = floor(size(cell_struct(1).(fd_names{1}),1)/bin_size); % trial length after binning
state_dim = setdiff([1:3],[time_dim, trial_dim]);
for ii = 1:length(cell_idx)
    i = cell_idx(ii);
    this_seq = [];
    num_trials = 0;
    for f = 1:numel(fd_names)
        try
            X = cell_struct(i).(fd_names{f});

            if ~isempty(fd_trial_idx) % only get specified trials
                if ~isempty(fd_trial_idx{f})
                    X = X(:,fd_trial_idx{f});
                end
            end
            this_num_trials = size(X,trial_dim);
            X = permute(X,[time_dim,trial_dim,state_dim]);
            if IF_MEDFILT
                for t = 1:this_num_trials
                    X(:,t) = medfilt1( X(:,t) ,3);
                end
            end
            
            % take average across specified frames if given
            if ~isempty(avg_frames)
                X = mean(X(avg_frames,:,:),1);
            end
            
            
            % dealing with more than one state
            XX = [];
            for tr = 1:this_num_trials
                XX = [XX;X(:,tr,:)];
            end
            XX = squeeze(XX);
            

            % binning
            if bin_size>1
                XX = cell2mat(arrayfun(@(x)mean(X((x-1)*bin_size+[1:bin_size],:,1)),1:num_bins,'UniformOutput', false)');
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
    this_num_trials =  size(cell_struct(1).(this_fd),trial_dim);
    if ~isempty(fd_trial_idx) % only get specified trials
        if ~isempty(fd_trial_idx{f})
            this_num_trials =numel(fd_trial_idx{f});
        end
            end
    catch
        this_num_trials = 0;
    end
    trial_idx.(this_fd) = trial_count+ [1:this_num_trials];
    trial_count = trial_count + this_num_trials;
end
end

