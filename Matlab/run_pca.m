function [ traj_struct, trial_idx,all_traj_struct] = run_pca( cell_struct,cell_idx_struct,opt,varargin )
train_struct = []; % default - compute svd
IF_USE_TRIAL_AVG = false; % use trial average to compute eigenvector 
IF_USE_WITHOLD_FRAMES = true; % only use frames in withold window to compute eigenvector
IF_USE_CORRECT_TRIALS = false; % only use correct trials to compute eigenvector
IF_USE_TIME_AVG = false; % take the average value within selected time window in a trial
IF_NORMALISE = false;


% opt.fd_names_avg is trial average traces
% opt.fd_names is individual trial traces
for v = 1:numel(varargin)
    if strcmpi(varargin{v},'train_struct')
        train_struct = varargin{v+1};
    end
    if strcmpi(varargin{v},'IF_USE_TRIAL_AVG')
        IF_USE_TRIAL_AVG = varargin{v+1};
    end
    if strcmpi(varargin{v},'IF_USE_CORRECT_TRIALS')
        IF_USE_CORRECT_TRIALS = varargin{v+1};
    end
    if strcmpi(varargin{v},'IF_USE_WITHOLD_FRAMES')
        IF_USE_WITHOLD_FRAMES = varargin{v+1};
    end
    if strcmpi(varargin{v},'IF_USE_TIME_AVG')
        IF_USE_TIME_AVG = varargin{v+1};
    end
    if strcmpi(varargin{v},'IF_NORMALISE')
        IF_NORMALISE = varargin{v+1};
    end
end
m_pca = opt.m; 
save_components = opt.save_components;
if ~IF_USE_WITHOLD_FRAMES
    trace_frames = 1:opt.trial_frames;
else
    trace_frames = opt.withold_frames_adj;
end
all_traj_struct = struct();

% get trial idx for cropping concatinated traces
trial_idx = struct();
trial_count = 0;
for f = 1:numel(opt.fd_names)
    this_fd = opt.fd_names{f}; 
    this_num_trials =  size(cell_struct(1).(this_fd),2);
    trial_idx.(this_fd) = trial_count+ [1:this_num_trials];
    trial_count = trial_count + this_num_trials;
end

% loop through cell_idx_fields (specified)
for idx = 1:numel(opt.idx_fields) % using different subsets of neuronsthis_idx_field = pa_opt.idx_fields{idx};
    this_idx_field = opt.idx_fields{idx};
    filt_cell_struct = cell_struct(cell_idx_struct.(this_idx_field));

    full_traces_in = []; %[trials, cells]
    traces_idx = struct();
    for f = 1:numel(opt.fd_names)
        this_traces = [];
        for c = 1:size(filt_cell_struct,2)
            this_trace = filt_cell_struct(c).(opt.fd_names{f})(:);
            this_traces = [this_traces this_trace];
        end
        temp_size = size(full_traces_in,1);
        full_traces_in = [full_traces_in; this_traces];
        traces_idx.(opt.fd_names{f}) = temp_size+1:size(full_traces_in,1);
    end
    full_traces_in = full_traces_in - mean(full_traces_in,2);
    svd_traces_in = full_traces_in; % use full trials for svd by default [time. cell]
    num_cells = size(svd_traces_in,2);
    
    if ~ IF_USE_CORRECT_TRIALS
        fd_names_svd = opt.fd_names;
    else
        fd_names_svd = opt.fd_names;
        fd_names_svd = fd_names_svd(contains(fd_names_svd,'_correct'));
    end
    
    
    if IF_USE_TIME_AVG  
        svd_traces_in = [];
        for f = 1:numel(opt.fd_names)
            this_traces = [];
            for c = 1:size(filt_cell_struct,2) % loop thro cells
                this_trace = filt_cell_struct(c).(opt.fd_names{f})(trace_frames,:);
                this_traces = [this_traces mean(this_trace,1)'];
            end
            svd_traces_in = [svd_traces_in; this_traces];
            
        end
        svd_traces_in = svd_traces_in - mean(svd_traces_in,2);
    end
    
    if IF_USE_TRIAL_AVG 
        avg_traces_in = [];
        for f = 1:numel(opt.fd_names_avg)
            this_traces = [];
            for c = 1:size(filt_cell_struct,2)
                this_trace = filt_cell_struct(c).(fd_names_svd)(:);
                this_concat = this_trace(trace_frames);
                this_traces = [this_traces this_concat];
            end
            temp_size = size(avg_traces_in,1);
            avg_traces_in = [avg_traces_in; this_traces];
            traces_idx.(opt.fd_names_avg{f}) = temp_size+1:size(avg_traces_in,1);
        end
        avg_traces_in = avg_traces_in - mean(avg_traces_in,2);
        svd_traces_in = avg_traces_in;
    end
    
    if IF_NORMALISE
        for c = 1:num_cells
        svd_traces_in(:,c) =  svd_traces_in(:,c)./max(svd_traces_in(:,c));
        end
    end

    if isempty(train_struct)
        if opt.UseGPU % takes longer than not using GPU.
            disp('Using GPU for svd...')
            tic
            gpu_trace = gpuArray(svd_traces_in');
            [U,S,V] = svd(gpu_trace);
            U = gather(single(U));
            S = gather(single(S));
            toc
        else
            tic
            [U,S,V] = svd(svd_traces_in'); % refer to matlab pca function
            toc
        end
        sigma_svd = diag(S);
        DOF = size(svd_traces_in,1)-1;
        latent_svd = sigma_svd.^2./DOF; % this is eigenvalues
        explained_svd = 100*latent_svd/sum(latent_svd);
        
    else
        U =  train_struct.U;
        explained_svd = [];
    end
    % principal component: eigenvector multiplied with sqrt of eigenvalue
    % P = US/sqrt(n-1); X = USV'; where U and V are orthogonal (UU' = I, VV' = I)
    save_components = 1:min([save_components(end),size(U,2)]); % PCs scaled by sqrt(n-1)
    m_pca = min([m_pca,save_components(end)]);
    
    F = full_traces_in*U(:,save_components);
    traj_struct.F = F;
    traj_struct.U = U;
    traj_struct.transmat = U(:,save_components);
    traj_struct.explained_svd = explained_svd;
    traj_struct.traces_in = full_traces_in;
    
    all_traj_struct.(this_idx_field) = traj_struct;
end



end
