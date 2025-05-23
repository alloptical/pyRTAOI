function [dpca_struct,traces_idx_struct,trial_idx_struct] = get_dpca_traj_brief(cell_struct,cell_idx,trial_fds,opt,varargin)
filt_cell_struct = cell_struct(cell_idx);
fig_save_path = opt.save_path;
decision_types = {'_correct','_incorrect'};
frame_range = 1:opt.trial_length;
IF_PLOT = true;
for v = 1:numel(varargin)
    if strcmpi(varargin{v},'IF_PLOT')
        IF_PLOT = varargin{v+1};
    end
    
     if strcmpi(varargin{v},'frame_range')
        frame_range = varargin{v+1};
    end
end

try
    save_name_ext= opt.save_name_ext;
catch
    save_name_ext = [];
end
N = size(filt_cell_struct,2);   % number of neurons
T = numel(frame_range);           % number of time points
S = 2;                          % number of stimuli
D = 2;                          % number of decisions
num_dPC = getOr(opt,'num_dPC',4);                   % number of PCs
num_trials = [];
for c = 1:numel(trial_fds)
    num_trials = [num_trials size(filt_cell_struct(1).(trial_fds{c}),2)];
end
maxTrialNum = max(num_trials);

firingRates = nan(N,S,D,T,maxTrialNum);
trialNum = nan(1,S,D);
for tex_count = 1:S
    for decision_count = 1:D
        this_fd = trial_fds(contains(trial_fds,decision_types{decision_count}));
        this_fd = this_fd{contains(this_fd,num2str(tex_count))};
        this_num_trials = size(filt_cell_struct(1).(this_fd),2);
        this_trial_traces = nan(N,T,maxTrialNum);
        for c = 1:N
            this_trial_traces(c,:,1:this_num_trials) = filt_cell_struct(c).(this_fd)(frame_range,:);
        end
        firingRates(:,tex_count,decision_count,:,:) = this_trial_traces; %data(neuron_count,frame_count,trial_count);
        trialNum(1,tex_count,decision_count) = this_num_trials;
    end
end
trialNum = repmat(trialNum,N,1);
firingRatesAverage = bsxfun(@times, nanmean(firingRates,5), size(firingRates,5)./trialNum);

% check consistency between trialNum and firingRates
for n = 1:size(firingRates,1)
    for s = 1:size(firingRates,2)
        for d = 1:size(firingRates,3)
            assert(isempty(find(isnan(firingRates(n,s,d,:,1:trialNum(n,s,d))), 1)), 'Something is wrong!')
        end
    end
end

combinedParams = {{1, [1 3]}, {2, [2 3]}, {3}, {[1 2], [1 2 3]}};
margNames = {'Stimulus', 'Decision', 'Condition-independent', 'S/D Interaction'};
margColours = [23 100 171; 187 20 25; 150 150 150; 114 97 171]/256;
tic

%  dPCA with regularization
ifSimultaneousRecording = 1;
optimalLambda = dpca_optimizeLambda(firingRatesAverage, firingRates, trialNum, ...
    'combinedParams', combinedParams, ...
    'simultaneous', ifSimultaneousRecording, ...
    'numRep', 10, ...  % increase this number to ~10 for better accuracy
    'filename', 'tmp_optimalLambdas.mat',...
    'numComps',4);

Cnoise = dpca_getNoiseCovariance(firingRatesAverage, ...
    firingRates, trialNum, 'simultaneous', ifSimultaneousRecording);

[W,V,whichMarg] = dpca(firingRatesAverage, num_dPC, ...
    'combinedParams', combinedParams, ...
    'lambda', optimalLambda, ...
    'Cnoise', Cnoise);
toc


close
%% dPCA decoding (separate components)
decodingClasses = {[(1:S)' (1:S)'], repmat([1:2], [S 1]), [], [(1:S)' (S+(1:S))']};
decodingClasses = {[1 1; 2 2], [1 2; 1 2], [], [1 2; 3 4]}; % test

try
    componentsSignif = [];
    if IF_PLOT
        % computing explained variance
        explVar = dpca_explainedVariance(firingRatesAverage, W, V, ...
            'combinedParams', combinedParams);
        
        accuracy = dpca_classificationAccuracy(firingRatesAverage, firingRates, trialNum, ...
            'lambda', optimalLambda, ...
            'combinedParams', combinedParams, ...
            'decodingClasses', decodingClasses, ...
            'simultaneous', ifSimultaneousRecording, ...
            'numRep', 10, ...        % increase to 100
            'filename', 'tmp_classification_accuracy.mat');
        
        accuracyShuffle = dpca_classificationShuffled(firingRates, trialNum, ...
            'lambda', optimalLambda, ...
            'combinedParams', combinedParams, ...
            'decodingClasses', decodingClasses, ...
            'simultaneous', ifSimultaneousRecording, ...
            'numRep', 10, ...        % increase to 100
            'numShuffles', 50, ...  % increase to 100 (takes a lot of time)
            'filename', 'tmp_classification_accuracy.mat');
        % try
        %     dpca_classificationPlot(accuracy, [], accuracyShuffle, [], decodingClasses,...
        %         'marginalizationNames', margNames,...
        %         'whichMarg', whichMarg,...
        %         'time', 1:opt.trial_length,...
        %         'timeEvents',opt.sta_gocue_frame)
        %
        %     suptitle([ opt.exp_name  this_idx_field ' ROIs' ' dPCA decoding accuracy'])
        %     fig=gcf;
        %     set(fig,'PaperPositionMode','auto');
        %     set(fig,'PaperOrientation','landscape');
        %     export_fig([fig_save_path filesep [opt.exp_name ' dPCA decoding accuracy ' this_idx_field ' ROIs'] '.png'])
        % end
        
        componentsSignif = dpca_signifComponents(accuracy, accuracyShuffle, whichMarg);
        dpca_plot(firingRatesAverage, W, V, @dpca_plot_default, ...
            'explainedVar', explVar, ...
            'marginalizationNames', margNames, ...
            'marginalizationColours', margColours, ...
            'whichMarg', whichMarg,                 ...
            'time', frame_range,                        ...
            'timeEvents',opt.sta_gocue_frame,               ...
            'timeMarginalization', 3,           ...
            'legendSubplot', 16,                ...
            'componentsSignif', componentsSignif);
        this_title = [ opt.exp_name  ' dPCA' save_name_ext];
        suptitle(strrep(this_title,'_', ' '))
        fig=gcf;
        set(fig,'PaperPositionMode','auto');
        set(fig,'PaperOrientation','landscape');
        export_fig([fig_save_path filesep this_title '.png'])
    end
    %% get dpca struct: trajectories and decoding mat
    % get the PCs that explain stimulus type
    PCidx.stim = find(whichMarg==1);
    PCidx.decision = find(whichMarg == 2);
    
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
    cell_mean =  mean(traces_in,1);
    traces_in = traces_in - cell_mean; % normalise to mean
    
    F = traces_in*W;
    traj_struct = struct();
    for f = 1:numel(trial_fds)
        this_trial_idx = traces_idx_struct.(trial_fds{f});
        this_traces = F(this_trial_idx,:);
        this_num_trials = size(this_traces,1)/opt.trial_length;
        this_traces_reshape = nan(this_num_trials,opt.trial_length,num_dPC); %[trial,frames,factor]
        for t = 1:this_num_trials
            for m = 1:num_dPC
                this_traces_reshape(t,:,m) = this_traces([1:opt.trial_length]+(t-1)*opt.trial_length,m)';
            end
        end       
            traj_struct.(trial_fds{f}) = this_traces_reshape;
    end
    % trial indices
    temp_num_trials = structfun(@(x)round(numel(x)./opt.trial_length),traces_idx_struct,'UniformOutput', false );
    trial_count = 0;
    for i = 1:numel(trial_fds)
        trial_idx_struct.(trial_fds{i}) = trial_count+[1:temp_num_trials.(trial_fds{i})];
        trial_count = trial_count+temp_num_trials.(trial_fds{i});
    end
    % save to structure
    this_dpca_struct.traj_struct = traj_struct;
    this_dpca_struct.W = W;
    this_dpca_struct.F = F;   
    this_dpca_struct.mean = cell_mean;
    this_dpca_struct.componentsSignif = componentsSignif;
    this_dpca_struct.whichMarg = whichMarg;
    this_dpca_struct.PCidx = PCidx;
    this_dpca_struct.U_stim = W(:,PCidx.stim);
    this_dpca_struct.U_decision = W(:,PCidx.decision);
    dpca_struct = this_dpca_struct;
    
catch e
    fprintf(1,'get dpca error! The message was:\n%s',e.message);
    return
end

end

