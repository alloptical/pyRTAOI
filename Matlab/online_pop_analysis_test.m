% plug-in for OnlineProcTexture
% parameters for running protocol PHOTO_WEIGH_SUM in pyRTAOI
% attempting to generate these params:
%     pop_params.weights;
%     pop_params.thresh;
%     pop_params.frames_enable_trigger;
%% load cell_struct
try
    [tex_file,tex_path] = uigetfile('*.mat','Select ProcTex data');
    disp(['Loaded ProcTex file :',fullfile(tex_path,tex_file)])
    load(fullfile(tex_path,tex_file));
catch
    warning('tex file not loaded')
end
%% init structures
cell_struct = tex_output_struct.cell_struct;
cell_idx_struct = tex_output_struct.cell_idx_struct;
cell_idx_struct.all = 1:size(cell_struct,2);
[trial_color] = online_tex_init_color();
%% ========================== SET PARAMETERS ==============================
opt = tex_output_struct.opt; % this will be shared across different analysis
opt.trial_color = trial_color;
opt.trigger_idx_fd = 'all'; % cell type used for analysis (a field in cell_idx_struct)
opt.trial_length= opt.sta_pre_frames+opt.sta_post_frames+1;
%% get fraction correct and incorrect
pc_correct.stim1 = size(cell_struct(1).('st_correct_stim_1'),2)/(size(cell_struct(1).('st_correct_stim_1'),2)+size(cell_struct(1).('st_incorrect_stim_1'),2));
pc_correct.stim2 = size(cell_struct(1).('st_correct_stim_2'),2)/(size(cell_struct(1).('st_correct_stim_2'),2)+size(cell_struct(1).('st_incorrect_stim_2'),2));
pc_incorrect.stim1 = 1-pc_correct.stim1;
pc_incorrect.stim2 = 1-pc_correct.stim2;
%% get epoch frames 
num_bs_frames = opt.sta_baseline_frames + 30;
num_stim_frames = 60; % number frames after baseline
num_decision_frames = 60;
epoch_frames.bs = 1:opt.sta_baseline_frames + 30;
epoch_frames.stim = epoch_frames.bs(end)+[1:num_stim_frames];
epoch_frames.decision = epoch_frames.stim(end)+[1:num_decision_frames];
%% HMM params (SEE HMMMAR wiki)
hmm_opt = struct();
hmm_opt.bin_size = 1; % binsize = 1 gave better result than 3; maybe because more data points
hmm_opt.fd_names = {'st_correct_stim_1','st_correct_stim_2','st_incorrect_stim_1','st_incorrect_stim_2'};
hmm_opt.num_states = 5;
hmm_opt.Fs = opt.frame_rate/hmm_opt.bin_size;            % sampling frequency
hmm_opt.K = 5;  	         % The number of states to infer
hmm_opt.order = 0; 	         % The lag used, this is only relevant when using MAR observations
hmm_opt.zeromean = 0; 	     % We do not want to model the mean, so zeromean is set on
hmm_opt.covtype = 'uniquefull';    % We want to model the full covariance matrix
hmm_opt.embeddedlags = 0; % 15 lags are used from -7 to 7
hmm_opt.pca = 0;   % The PCA dimensionality reduction is 2 times the number of ROIs because the covar matrix bins x channels is high
hmm_opt.initrep = 1; % to make it quicker - leave by default otherwise
hmm_opt.initcyc = 1; % to make it quicker - leave by default otherwise 
hmm_opt.cyc = 200; % to make it quicker - leave by default otherwise
hmm_opt.verbose = 1; 
hmm_opt.standardise = 0;
hmm_opt.go_cue_bin = round(opt.sta_gocue_frame/hmm_opt.bin_size);
hmm_opt.idx_field = opt.trigger_idx_fd ; % idx field name in cell_idx_struct 
hmm_opt.max_num_run = 10; % maximum number of training cycle before states are related with trial structure 
%% PCA params
pca_opt.UseGPU = false;
pca_opt.m = 3;
pca_opt.Fs = opt.frame_rate;
pca_opt.save_components = pca_opt.m;
pca_opt.trial_frames = opt.sta_pre_frames+opt.sta_post_frames+1;
pca_opt.withold_frames_adj = epoch_frames.stim;
pca_opt.idx_fields = {opt.trigger_idx_fd};
pca_opt.go_cue_bin = opt.sta_pre_frames;
pca_opt.fd_names = {'st_correct_stim_1','st_correct_stim_2','st_incorrect_stim_1','st_incorrect_stim_2'};
pca_opt.trial_length = opt.trial_length;
%% FA params
fa_opt = struct();
fa_opt.m = 5;
fa_opt.bin_size = 1;
fa_opt.go_cue_bin = round(opt.sta_gocue_frame/fa_opt.bin_size);
fa_opt.Fs = opt.frame_rate/fa_opt.bin_size;            % sampling frequency
fa_opt.trial_length = floor(opt.trial_length/fa_opt.bin_size);
fa_opt.save_components = pca_opt.m;
fa_opt.trial_frames = opt.sta_pre_frames+opt.sta_post_frames+1;
fa_opt.withold_frames_adj = epoch_frames.stim;
fa_opt.idx_fields = {opt.trigger_idx_fd};
fa_opt.go_cue_bin = opt.sta_pre_frames;
fa_opt.fd_names = {'st_correct_stim_1','st_correct_stim_2','st_incorrect_stim_1','st_incorrect_stim_2'};

%% Coding derectio params
cd_opt = fa_opt;
cd_opt.idx_fields = {'all'};
cd_opt.frames_to_avg = epoch_frames.decision;

%% ================  GET POPULATION TRAJECTORIES ==========================
%% PCA
disp('Rnning PCA ...')
plot_fields = pca_opt.fd_names;
[pca_struct,pca_trial_idx ] = run_pca( cell_struct,cell_idx_struct,pca_opt);
[pca_traj_struct,num_pcs] = get_pop_vectors(pca_struct.F,pca_opt.trial_frames,pca_trial_idx);
plot_pop_vectors(pca_traj_struct,plot_fields,pca_opt.m,pca_opt,'plot_ylabel','PC level');
disp('...Done')

%% FA
disp('Rnning factor analysis...')
[fa_traces_in,fa_trial_idx] = get_input_seq(cell_struct,cell_idx_struct.(fa_opt.idx_fields{1}),fa_opt.fd_names,fa_opt.bin_size);%%
fa_struct = struct();
fa_struct.mean = mean(fa_traces_in);
fa_struct.std = std(fa_traces_in);
[fa_struct.lambda,fa_struct.psi,fa_struct.T,fa_struct.stats,fa_struct.F] = factoran(fa_traces_in,fa_opt.m,'Xtype','data','Maxit',1000);
sqrtPsi = sqrt(fa_struct.psi); % get transition matrix (multiply to zscored data)
invsqrtPsi = diag(1 ./ sqrtPsi);
fa_struct.transmat = invsqrtPsi/(fa_struct.lambda'*invsqrtPsi);
[fa_traj_struct,num_fs] = get_pop_vectors(fa_struct.F,fa_opt.trial_length,fa_trial_idx);
plot_pop_vectors(fa_traj_struct,plot_fields,fa_opt.m,fa_opt,'plot_ylabel','Factor level');
disp('...Done')

%% CD
disp('Getting coding directions...')
[cd_struct,this_cd,this_db,cd_struct_framewise,this_framewise_cd] = run_cd(cell_struct(cell_idx_struct.(cd_opt.idx_fields{1})),cd_opt);
cd_struct.transmat = this_cd;
plot_pop_vectors(cd_struct,plot_fields,1,cd_opt,'plot_ylabel','Coding direction');
disp('...Done')

%% HMM
% https://github.com/OHBA-analysis/HMM-MAR/wiki/
% a MAR model aims to predict the multichannel signal value at each time point as a linear combination of signal values at previous time points.
% organise data for hmmmar (each row is a sequence)
disp('Rnning HMM analysis...')
epoch_bins = structfun(@(x)unique(ceil(x./hmm_opt.bin_size)),epoch_frames,'UniformOutput',false);
num_bins = floor(size(cell_struct(1).(hmm_opt.fd_names{1}),1)/hmm_opt.bin_size); % trial length after binning
[seq,trial_idx,num_trials] = get_input_seq(cell_struct,cell_idx_struct.('all'),hmm_opt.fd_names,hmm_opt.bin_size);%%
% run unsupervised hmm-mar 
% train using forward-backward algorithm 
T = num_bins * ones(num_trials,1);
datain = seq;
epoch_window = epoch_bins.decision;
ref_fields = {'st_correct_stim_1','st_correct_stim_2'};
plot_fields = hmm_opt.fd_names;
run_count = 0;
FLAG_MATCHED_STATES = 0;
while (~(all(FLAG_MATCHED_STATES>0))) && (run_count<hmm_opt.max_num_run)
    run_count = run_count+1;
    [hmm,Gamma_hmm, ~, vpath] = hmmmar(datain,T,hmm_opt);
    [hmm_struct,num_states] = get_pop_vectors(Gamma_hmm,num_bins,trial_idx);
    % repeat if no state was found correlated with trial structure
    [state_idx_struct,corr_struct,FLAG_MATCHED_STATES] = get_state_idx(hmm_struct,epoch_window,ref_fields,num_states);
    disp(['Number HMM run: ' num2str(run_count)])
end
hmm_opt.m = num_states;
disp('...Done')
plot_pop_vectors(hmm_struct,plot_fields,num_states,hmm_opt);

%% show Viterbi path
vpath_struct = get_pop_vectors(vpath,num_bins,trial_idx);
plot_pop_vectors(vpath_struct,plot_fields,1,hmm_opt,'ylimit',[0, hmm_opt.K],'plot_ylabel','Viterbi path');

%% number of state transitions before go-cue
fds = hmm_opt.fd_names;
frames_of_interest = [-30:1:0]+ opt.sta_pre_frames;
[trans_struct] = plot_num_trans(vpath_struct,frames_of_interest,fds,trial_color);
%% get steady state frames from viterbi path
search_window = epoch_bins.stim;
windows = struct();
fds = hmm_opt.fd_names;
for f = 1:numel(fds)
    fd = fds{f};
    this_traces = vpath_struct.(fd);
    this_windows = zeros(size(this_traces));
    for t = 1:size(this_traces,1)
        X = this_traces(t,:);
        idx = find(X == state_idx_struct.(strrep(fd,'_incorrect_','_correct_')));
        idx = intersect(idx, search_window);
        this_windows(t,idx) = 1;
    end
    windows.(fd) = this_windows;
    figure
    imagesc(this_windows)
    title(fd)
end
%% ======================= GET DECODERS =========================
% select which population vector to use
opt.decod_type = 'fa'; % fa, cd, hmm or pca
opt.IF_PROJ_TO_STIM = 0;
opt.Nstd = 2;
switch opt.decod_type
    case 'pca'
        this_traj_struct = pca_traj_struct;
        this_struct = pca_struct;
        this_opt = pca_opt;
        opt.IF_NORMALISE = 0;
    case 'fa'
        this_traj_struct = fa_traj_struct;
        this_struct = fa_struct;
        this_opt = fa_opt;
        opt.IF_NORMALISE = 1; % need to take zscore if using factor analysis

    case 'hmm'
        this_traj_struct = hmm_struct;
        this_struct = hmm;
        this_opt = hmm_opt;
        opt.IF_NORMALISE = 0;
    case 'cd'
        this_traj_struct = cd_struct;
        this_struct = cd_struct;
        this_opt = cd_opt;  
        opt.IF_NORMALISE = 0;
end
this_opt.frames_to_avg = [-30:1:-1]+this_opt.go_cue_bin;
this_opt.gocue_frame_adj = this_opt.go_cue_bin;
this_opt.IF_FRAMEWISE = true;
this_opt.Nstd = opt.Nstd;
this_opt.min_frames = 5;
this_opt.disc_frames_fd = 'disc_frames';
this_opt.disc_frame_fd = 'disc_frame';
this_opt.trial_color = opt.trial_color;
%% TRAIN CLASSIFIER stim1 vs stim2
if (~strcmp( opt.decod_type,'cd')) && opt.IF_PROJ_TO_STIM
    [stim_classifier] = get_stimtype_classifier(this_traj_struct,windows,this_opt);
    disp(['Stim type classifier accuracy ' num2str(stim_classifier.accur)])
    proj_fds = this_opt.fd_names;
    weights = stim_classifier.B(2:end);
    [proj_struct] = get_projections(this_traj_struct,weights,proj_fds);
    
else
    proj_struct = this_traj_struct;
end
% Project trial data onto stim type classifier
plot_pop_vectors(proj_struct,plot_fields,1,this_opt)

%% TRAIN CLASSIFIER correct vs incorrect
decod_struct = struct();
decod_struct =  get_trialoutcome_classifier( decod_struct,proj_struct, this_opt,'IF_FRAMEWISE',this_opt.IF_FRAMEWISE,'fd_names',this_opt.fd_names);
if ~opt.IF_PROJ_TO_STIM % project to outcome classifiers
    [proj_struct] = get_projections(proj_struct,decod_struct.B.stim1(2:end),{'st_correct_stim_1','st_incorrect_stim_1'},'proj_struct',proj_struct);
    [proj_struct] = get_projections(proj_struct,decod_struct.B.stim2(2:end),{'st_correct_stim_2','st_incorrect_stim_2'},'proj_struct',proj_struct);
end
[ decod_struct ] =  get_disc_time( proj_struct, decod_struct,this_opt,'IF_FRAMEWISE',this_opt.IF_FRAMEWISE);
figure('name','trial outcome decoder')
plot_decoder_performance(decod_struct,proj_struct,[],trial_color,this_opt);

%% ==================== GENERATE OUTPUT STRUCTURE =========================
% in pyRTAOI photo_frames = [wait_frames monitor_frames] - CHECK THE
% SETTINGS!
condition_type = 2; % change this
opt.pop_opt = this_opt;
pop_params = struct();
pop_params.condition_type = condition_type;
if ~strcmp( opt.decod_type,'cd')&& opt.IF_PROJ_TO_STIM
    pop_params.weights = this_struct.transmat* stim_classifier.B(2:end);
elseif ~opt.IF_PROJ_TO_STIM
    pop_params.weights = this_struct.transmat*decod_struct.B.(['stim' num2str(condition_type)])(2:end);
elseif strcmp( opt.decod_type,'cd') 
    pop_params.weights = this_struct.transmat;
end
pop_params.thresh = decod_struct.thresh_fix.(['stim' num2str(condition_type)]);
if(opt.IF_NORMALISE )
    [pop_params.weights,pop_params.thresh] = get_norm_weights(pop_params.weights,pop_params.thresh,this_struct.mean,this_struct.std);
end

pop_params.frames_enable_trigger = [0:1:5] + decod_struct.(['stim' num2str(condition_type)]).(this_opt.disc_frame_fd);
[output] = generate_cell_idx_file(cell_struct,cell_idx_struct,pop_params,opt);


%% ============ test decoder performance on testset ======================= 
this_opt = output.pop_opt;
test_decod_struct = struct();
photostim_frames = output.trigger_frames
thresh_struct.stim1 = output.trigger_thresh;
thresh_struct.stim2 = output.trigger_thresh;

[test_traces_in,test_trial_idx] = get_input_seq(cell_struct,output.trigger_idx,this_opt.fd_names,this_opt.bin_size);%%
test_proj_struct = get_pop_vectors(test_traces_in,this_opt.trial_length,test_trial_idx);
[test_proj_struct] = get_projections(test_proj_struct,output.trigger_weights,this_opt.fd_names);
[ test_decod_struct ] =  get_disc_time( test_proj_struct, test_decod_struct,this_opt,'threshold',thresh_struct,...
    'IF_FRAMEWISE',this_opt.IF_FRAMEWISE,'IF_USE_CROSSVAL_RESULT',false);
[test_decod_struct] = detect_error_trials(test_decod_struct,test_proj_struct,thresh_struct,[],photostim_frames,this_opt);
figure('name','trial outcome decoder test')
plot_decoder_performance(test_decod_struct,test_proj_struct,[],trial_color,this_opt);













%% decode states from new data using trained hmm
% using same data as training set just to test
hmm.algorithm = 'forward'; % test using forward algorithm 
[Gamma,Xi] = hmmdecode(datain,T,hmm,0,[],0);
[viterbipath] = hmmdecode(datain,T,hmm,1,[],0);
onsets = getStateOnsets(viterbipath,T,hmm_opt.Fs,hmm.K);
plot_pop_vectors(Gamma,num_bins,trial_idx,hmm_opt)



%% show onset time - looks same cross conditions, should plot transition times
num_states = size(onsets,2);
onset_times = struct();
fds = fields(trial_idx);
for f = 1:numel(fds)
    trial_num = trial_idx.(fds{f});
    this_onset = [];
    for s = 1:size(onsets,2)
        for i = 1:length(trial_num)
            this_onset = [this_onset;onsets{i,s}];
        end
        onset_times(s).(fds{f}) = this_onset;
    end
end
figure('name','state onset times')
state_colors = getOr(opt,'state_colors',brewermap(num_states,'Set1'));
condi_colors = cell2mat(cellfun(@(f)trial_color.(f),fds,'UniformOutput',false));
for s = 1:num_states
    subplot(num_states,1,s)
    scatter_cmp_conditions(onset_times(s),[],...
        1,condi_colors,'connect_scatter',0,'BriefXlabel',1,'ShowMeanInXlabel',1,'add_jitter',1);
    title(['State ' num2str(s)])
    ylabel(['Onset time(s)'])
end












% ==================== STUFF TRIED BUT DID NOT WORK =======================
%% initialise transition matrix by performance - NOT USED FOR HMM-MAR
% states: baseline; stim1; stim2; decision1, decision2
num_bins_bs = num_bs_frames/hmm_opt.bin_size;
num_bins_stim = num_stim_frames/hmm_opt.bin_size;
num_bins_decision = num_decision_frames/hmm_opt.bin_size;

initTR = [1-1/num_bins_bs, 1/num_bins_bs*0.5,1/num_bins_bs*0.5, 0, 0;...
          0,1-1/num_bins_stim,0,1/num_bins_stim*pc_correct.stim1,1/num_bins_stim*pc_incorrect.stim1;...
          0,1-1/num_bins_stim,0,1/num_bins_stim*pc_incorrect.stim2,1/num_bins_stim*pc_correct.stim2;...
          1/num_bins_decision,0,0,1-1/num_bins_decision,0;...
          1/num_bins_decision,0,0,0,1-1/num_bins_decision];
%% initialise emission matrix by average activity during epochs - NOT USED FOR HMM-MAR

cell_avg.bs = cell2mat(arrayfun(@(x)mean([mean(x.('st_correct_stim_1')(epoch_frames.bs,:),1),mean(x.('st_correct_stim_2')(epoch_frames.bs,:),1)]),cell_struct,'UniformOutput' , false));
cell_avg.stim1 = cell2mat(arrayfun(@(x)mean([mean(x.('st_correct_stim_1')(epoch_frames.stim,:),1),mean(x.('st_incorrect_stim_1')(epoch_frames.stim,:),1)]),cell_struct,'UniformOutput' , false));
cell_avg.stim2 = cell2mat(arrayfun(@(x)mean([mean(x.('st_correct_stim_2')(epoch_frames.stim,:),1),mean(x.('st_incorrect_stim_2')(epoch_frames.stim,:),1)]),cell_struct,'UniformOutput' , false));
cell_avg.decision1 = cell2mat(arrayfun(@(x)mean([mean(x.('st_correct_stim_1')(epoch_frames.decision,:),1),mean(x.('st_incorrect_stim_2')(epoch_frames.decision,:),1)]),cell_struct,'UniformOutput' , false));
cell_avg.decision2 = cell2mat(arrayfun(@(x)mean([mean(x.('st_correct_stim_2')(epoch_frames.decision,:),1),mean(x.('st_incorrect_stim_1')(epoch_frames.decision,:),1)]),cell_struct,'UniformOutput' , false));

initE = [cell_avg.bs; cell_avg.stim1;cell_avg.stim2;cell_avg.decision1;cell_avg.decision2];


%% run suppervised hmm-mar - dosen't work well, stick with unsuppervised
% decoding weights are same across trials; but decoding window is different
% make stimulus feature matrix
hmm_opt.covtype = 'uniquediag'; % dosen't work with full cov matrix
hmm_opt.K = 4;
hmm_opt.order = 1;
hmm_opt.fd_names_sup = {'st_correct_stim_1','st_correct_stim_2'}; % train only with correct trials
num_sup_trials = length( cell2mat(cellfun(@(x)trial_idx.(x),hmm_opt.fd_names_sup,'UniformOutput',false))); 
Tsup = (num_bins-1) * ones(num_sup_trials,1);
[concat_seq,trial_idx_sup] = get_input_seq(cell_struct,cell_idx_struct.('tex'),hmm_opt.fd_names_sup,hmm_opt.bin_size);%%

% three states corresponding to stim1, stim2 and baseline - dosent work,
% taking into account timecourse within trials
Yin = zeros(size(concat_seq,1),hmm.K);
stim_bin_idx = unique(round([epoch_frames.stim epoch_frames.decision]./hmm_opt.bin_size));
bs_bin_idx = unique(round( epoch_frames.bs./hmm_opt.bin_size));
for f = 1:numel(hmm_opt.fd_names_sup)
    this_fd = hmm_opt.fd_names_sup{f};
    this_trial_idx = trial_idx.(this_fd);
    for i = this_trial_idx
        Yin((i-1)*(num_bins-1)+1+stim_bin_idx,f ) = 1;
        Yin((i-1)*(num_bins-1)+1+bs_bin_idx,3) = 1;
    end
end

% mark trial types (stim 1 or stim2) without timecourse within trials
Yin = zeros(size(Tsup,1),1);
% specifying stimulus types
for f = 1:numel(hmm_opt.fd_names_sup)
    this_fd = hmm_opt.fd_names_sup{f};
    this_trial_idx = trial_idx.(this_fd);
    Yin(this_trial_idx) = (-1)^f;
end
[tuda_stim1_vs_stim2,Gamma_tuda_stim1_vs_stim2] = tudatrain(concat_seq,Yin,Tsup,hmm_opt);
tuda_struct = plot_pop_vectors(Gamma_tuda_stim1_vs_stim2,num_bins-hmm_opt.order,trial_idx_sup,hmm_opt);

%%
[data_test,trial_idx_test] = get_input_seq(cell_struct,cell_idx_struct.('tex'),hmm_opt.fd_names,hmm_opt.bin_size);%%
[Gamma,Xi] = hmmdecode(concat_seq,Tsup,tuda_stim1_vs_stim2,0,[],0);
[viterbipath] = hmmdecode(data_test,T,tuda_stim1_vs_stim2,1,[],0);
%% outcome decoder using  corret and incorrect trials (without projecting to stim decoder)
% project to outcome decoders directly
decod_direct= struct();
decod_direct =  get_trialoutcome_classifier( decod_direct,this_traj_struct, this_opt,'IF_FRAMEWISE',this_opt.IF_FRAMEWISE,'fd_names',this_opt.fd_names);
[proj_test] = get_projections(this_traj_struct,decod_direct.B.stim1(2:end),{'st_correct_stim_1','st_incorrect_stim_1'});
[proj_test] = get_projections(this_traj_struct,decod_direct.B.stim2(2:end),{'st_correct_stim_2','st_incorrect_stim_2'},'proj_struct',proj_test);
[ decod_direct ] =  get_disc_time( proj_test, decod_direct,this_opt,'IF_FRAMEWISE',this_opt.IF_FRAMEWISE,'IF_USE_CROSSVAL_RESULT',true);

figure('name','trial outcome decoder')
plot_decoder_performance(decod_direct,proj_test,[],trial_color,this_opt);
