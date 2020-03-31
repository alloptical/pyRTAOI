% Posthoc analysis on go-nogo pyrtaoi recordings
% 1. reproduce online plots
% 2. train decoders using trace in cell_struct rather than raw online trace
% 3. trial alternative ways of decoding choice
% adapted from OnlineProcTexture
% ZZ 2020

%% add path 
clear all
close all
clc

% ZZ PC
% addpath(genpath('C:\Users\User\Desktop\pyRTAOI-rig\Matlab'));
% cd('C:\Users\User\Desktop\pyRTAOI-rig\Matlab');

% BRUKER1
% matlab_set_paths_zz

IF_RECOMPUTE_STA = false;
IF_GO_NOGO = true;
opt = struct();
opt.IF_GO_NOGO = true;
%% load online saved data struct (ProxTex...)
try
%% load online proc data
[online_file,online_path] = uigetfile('*.mat','Select texture caiman data');
disp(['Loaded file :',fullfile(online_path,online_file)])
online_data = load(fullfile(online_path,online_file)); 
online_data = online_data.tex_output_struct;
caiman_path = online_path;

%% parse data
opt = online_data.opt;
trial_color = opt.trial_color;
cell_idx_struct = online_data.cell_idx_struct;
cell_struct = online_data.cell_struct;
pop_params = online_data.pop_params;
IF_ONLINE_DATA_LOADED = true;
catch
    IF_RECOMPUTE_STA = true;
    IF_ONLINE_DATA_LOADED = false;
    disp('no online struct loaded')
end
%% load caiman data (if no online saved data is loaded or if want to recompute sta)
% recompute sta using zscored full trace - dpca and stim decoder weights looks the same as unzscored
if IF_RECOMPUTE_STA||~IF_ONLINE_DATA_LOADED
    %% load CAIMAN data
    if IF_ONLINE_DATA_LOADED %get file name from online data
        [~,caiman_file] = fileparts(online_data.input_caiman_file);
        caiman_file = [online_path filesep caiman_file '.mat'];
        caiman_data = load(caiman_file);
        disp(['Loaded file :',caiman_file])
    else % select caiman file
        [caiman_file,caiman_path] = uigetfile('*.mat','Select caiman data');
        caiman_data = load(fullfile(caiman_path,caiman_file));
        disp(['Loaded file :',fullfile(caiman_path,caiman_file)])
        
    end
    %% load pybehav data
    try
        [pb_file,pb_path] = uigetfile('*.mat','Select pybehavior data',caiman_path);
        disp(['Loaded file :',fullfile(pb_path,pb_file)])
        behavior_data =  load(fullfile(pb_path,pb_file)) ;
        trials = make_trials_struct(behavior_data);
        FLAG_PYBEHAV_LOADED = true;
    catch
        FLAG_PYBEHAV_LOADED = false;
    end
    %% load config file: baseline output struct (file named as 'OutputParams')
    [baseline_file,baseline_path] = uigetfile('*.mat','Select baseline OutputParams',caiman_path);
    baseline_output = load(fullfile(baseline_path,baseline_file));
    baseline_output = baseline_output.output;
    disp(['Loaded file :',fullfile(baseline_path,baseline_file)])

    %% initialise parameters
    if ~IF_ONLINE_DATA_LOADED
        [opt] = init_opt_posthoc(opt);
        opt.discard_trials_after = [];
        opt.discard_trials_before = [];
        opt.exp_name = ['Posthoc_' baseline_output.exp_name];
        trial_color = online_tex_init_color();
    end
    %% sort trial types
    [trials,trial_indices,sens_stim_frames] = sort_trial_types(trials,caiman_data,FLAG_PYBEHAV_LOADED,IF_GO_NOGO,opt);
    [dp,hr,fa] = get_gonogo_baseline_performance(trial_indices);
    %% make a cell_idx_struct and pop_params based on output file
    if ~IF_ONLINE_DATA_LOADED
        cell_idx_struct.('tex') =  baseline_output.trigger_idx;
        cell_idx_struct.tex1 = baseline_output.target_ensembles{1};
        cell_idx_struct.tex2 = baseline_output.target_ensembles{2};
        cell_idx_struct.all = 1:num_cells;
        opt.target_idx_fd = {'tex1','tex2'};
        opt.trigger_idx_fd = 'tex';
        opt.method='dpca';
        pop_params.weights = baseline_output.trigger_weights;
    end
    %% get cnm and cell structures
    [cnm_struct,cell_struct] = get_cell_struct(caiman_data,sens_stim_frames,opt);
    %% sort STAs for different trial types
    trial_types = fields(trial_indices);
    raw_cell_struct = struct();
    for i = 1:numel(trial_types)
        this_fd = trial_types{i};
        this_idx = trial_indices.(this_fd);
        this_idx = this_idx(this_idx<num_trials);
        if ~isempty(this_idx)
            for c = 1:num_cells
                cell_struct(c).(this_fd) = cell_struct(c).sta_traces( this_idx,:)';
                raw_cell_struct(c).(this_fd) = cell_struct(c).raw_sta_traces( this_idx,:)';
            end
        else
            for c = 1:num_cells
                cell_struct(c).(this_fd) = [];
                raw_cell_struct(c).(this_fd) = [];
            end
        end
        
    end
    disp('sorted cell_struct sta_traces')
    %% add some fields to cell_struct
    if ~IF_ONLINE_DATA_LOADED
        for i = 1:num_cells
            if ~isempty(intersect(i,baseline_output.trigger_idx))
                cell_struct(i).is_tuned = 1;
            else
                cell_struct(i).is_tuned = 0;
            end
        end
    end

end
%% config save path
save_path = [caiman_path filesep 'analysis_files'];
if ~exist(save_path, 'dir')
    mkdir(save_path)
end

fig_save_path = [caiman_path filesep 'figures'];
if ~exist(fig_save_path, 'dir')
mkdir(fig_save_path)
end

opt.save_path = fig_save_path;
opt.exp_name = strrep(online_file,'.mat','');
opt.save_name_ext = 'posthoc_';
%% normalise to baseline 
fd_names = {'stim_1_correct','stim_1_incorrect',...
            'stim_2_correct','stim_2_incorrect',...
};
num_cells = size(cell_struct,2);
for i = 1:num_cells
    for f = 1:numel(fd_names)
        cell_struct(i).(fd_names{f}) =  cell_struct(i).(fd_names{f}) - mean(cell_struct(i).(fd_names{f})(1:opt.sta_baseline_frames,:),1);
    end
end
disp('normalised sta')
%% Plot STAs traces for certain cell types
plot_cell_type = 'tex'; % trigger cells
peak_frame_range = opt.sta_peak_search_range;
avg_frame_range = opt.sta_avg_frames;

plot_cell_idx = cell_idx_struct.(plot_cell_type);
% rank plot cells by weights
[~,temp_idx] = sort(pop_params.weights);
plot_cell_idx = plot_cell_idx(temp_idx);
plot_cell_w = pop_params.weights(temp_idx);

target_cell_idx = unique(cell2mat(cellfun(@(x)cell_idx_struct.(x),opt.target_idx_fd,'UniformOutput',false)));
num_plot_cols = 8;
num_plot_rows = ceil(numel(plot_cell_idx)/num_plot_cols);
num_stim_type = 2;

% raw traces
figure('name','tex cell sta traces','units','normalized','outerposition',[0 0 1 1])
for ii = 1:numel(plot_cell_idx)
    i = plot_cell_idx(ii);
    subtightplot(num_plot_rows,num_plot_cols,ii)
    hold on
    % correct trials
    for j = 1:num_stim_type
        plot(cell_struct(i).(['stim_' num2str(j) '_correct']),'color',trial_color.(['correct_stim' num2str(j)]),'linewidth',1)
        plot(cell_struct(i).(['stim_' num2str(j) '_incorrect']),'color',trial_color.(['incorrect_stim' num2str(j)]),'linewidth',1)
    end
    
    xlim([0 opt.sta_pre_frames+opt.sta_post_frames+1])
    axis square
    
    % mark targets
    if ~isempty(intersect(target_cell_idx,i))
       text(0,0,'TARGET','units','normalized','color',trial_color.photo,'Horizontalalignment','left','VerticalAlignment','bottom')
       this_tex_color = trial_color.photo;
    else
        this_tex_color = 'black';
    end

    text(1,1,['Cell ' num2str(i) ' (ROI ' num2str(cell_struct(i).cnm_idx) ')'],'units','normalized','color',this_tex_color,'Horizontalalignment','right','VerticalAlignment','top')
    text(1,0,['w = ' num2str(plot_cell_w(ii),'%0.2f')],'units','normalized','color',this_tex_color,'Horizontalalignment','right','VerticalAlignment','bottom')


    box off
    
    % mark tuned cells
    if( cell_struct(i).is_tuned)
        this_pref = 0;
        
        if ~isempty(intersect(cell_idx_struct.tex1,i))
            this_pref = 1;
        elseif  ~isempty(intersect(cell_idx_struct.tex2,i))
            this_pref = 2;
        end
        if this_pref>0
            box on
            this_color = trial_color.(['correct_stim' num2str(this_pref)]);
            set(gca,'XColor',this_color,'YColor',this_color,'linewidth',3)
        end
    end
    
    
    % mark time
    plot([1,1].*opt.sta_gocue_frame, ylim,'color','black','linestyle',':')
    if ~ opt.flag_use_peak
        plot([avg_frame_range(1),avg_frame_range(end)],[0,0],'color',[.5 .5 .5],'linewidth',2)
    else
        plot([peak_frame_range(1),peak_frame_range(end)],[0,0],'color',[.5 .5 .5],'linewidth',2)
    end
    
    % show auc
    try
    text(0.05,.8,['trial auc ' num2str(cell_struct(i).correct_stimAUC,'%10.2f')],'units','normalized', 'horizontalalignment','left', 'color','black')
    text(0.05,.7,['zscore auc '  num2str(cell_struct(i).correct_stimAUC_zscore,'%10.2f') ],'units','normalized', 'horizontalalignment','left', 'color','black')
    end
    
    % modify x ticks
    xaxisvalues = [0:30:opt.sta_pre_frames+opt.sta_post_frames];
    xticks(xaxisvalues)
    xticklabels(arrayfun(@(x){num2str(x)},(xaxisvalues-opt.sta_trialon_frame)./opt.frame_rate))
    
%     ylim([-5 20])
end
% export_fig([fig_save_path filesep plot_cell_type 'STATrace_' strrep(caiman_file,'.mat','.png')])

% shaded error bar
figure('name','shaded tex cell sta traces','units','normalized','outerposition',[0 0 1 1])
x_ticks =[0:1:opt.trial_length-1];

for  ii = 1:numel(plot_cell_idx)
    i = plot_cell_idx(ii);
    subtightplot(num_plot_rows,num_plot_cols,ii)
    hold on
    for j = 1:num_stim_type
        try
        this_traces = cell_struct(i).(['stim_' num2str(j) '_correct'])';
        shadedErrorBar(x_ticks,mean(this_traces,1),...
            std(this_traces,[],1),{'color',trial_color.(['correct_stim' num2str(j)]),'linewidth',2},0.1);
        end
        try
        this_traces = cell_struct(i).(['stim_' num2str(j) '_incorrect'])';
        shadedErrorBar(x_ticks,mean(this_traces,1),...
            std(this_traces,[],1),{'color',trial_color.(['incorrect_stim' num2str(j)]),'linewidth',2},0.1);
        end
    end
    
    xlim([0 opt.sta_pre_frames+opt.sta_post_frames+1])
    axis square
    
    % mark targets
    if ~isempty(intersect(target_cell_idx,i))
       text(0,0,'TARGET','units','normalized','color',trial_color.photo,'Horizontalalignment','left','VerticalAlignment','bottom')
       this_tex_color = trial_color.photo;
    else
        this_tex_color = 'black';
    end

    text(1,1,['Cell ' num2str(i) ' (ROI ' num2str(cell_struct(i).cnm_idx) ')'],'units','normalized','color',this_tex_color,'Horizontalalignment','right','VerticalAlignment','top')
    text(1,0,['w = ' num2str(plot_cell_w(ii),'%0.2f')],'units','normalized','color',this_tex_color,'Horizontalalignment','right','VerticalAlignment','bottom')


    box off  
    % mark tuned cells
    if( cell_struct(i).is_tuned)
        this_pref = 0;
        
        if ~isempty(intersect(cell_idx_struct.tex1,i))
            this_pref = 1;
        elseif  ~isempty(intersect(cell_idx_struct.tex2,i))
            this_pref = 2;
        end
        if this_pref>0
            box on
            this_color = trial_color.(['correct_stim' num2str(this_pref)]);
            set(gca,'XColor',this_color,'YColor',this_color,'linewidth',3)
        end
    end
    

    % mark time
    plot([1,1].*opt.sta_gocue_frame, ylim,'color','black','linestyle',':','linewidth',2)
    if ~ opt.flag_use_peak
        plot([avg_frame_range(1),avg_frame_range(end)],[0,0],'color',[.5 .5 .5],'linewidth',2)
    else
        plot([peak_frame_range(1),peak_frame_range(end)],[0,0],'color',[.5 .5 .5],'linewidth',2)
    end
    
    % show auc
    try
    text(0.05,.8,['trial auc ' num2str(cell_struct(i).correct_stimAUC,'%10.2f')],'units','normalized', 'horizontalalignment','left', 'color','black')
    text(0.05,.7,['zscore auc '  num2str(cell_struct(i).correct_stimAUC_zscore,'%10.2f') ],'units','normalized', 'horizontalalignment','left', 'color','black')
    end
    % modify x ticks
    xaxisvalues = [0:30:opt.sta_pre_frames+opt.sta_post_frames];
    xticks(xaxisvalues)
    xticklabels(arrayfun(@(x){num2str(x)},(xaxisvalues-opt.sta_trialon_frame)./opt.frame_rate))
    
end
% export_fig([fig_save_path filesep plot_cell_type 'ShadedSTATrace_' strrep(caiman_file,'.mat','.png')])

%% Plot trial average as imagesc - photo response
plot_photo_types = {'_correct','_incorrect'};
plot_condition_types = [1,2];
cnn_thresh = min(min(cell2mat({cell_struct(cell_idx_struct.(opt.trigger_idx_fd)).cnn_prediction})),0.1);
plot_cell_idx = find(cell2mat({cell_struct(:).cnn_prediction})>=cnn_thresh);
photo_ensembles{1} = cell_idx_struct.tex1;
photo_ensembles{2} = cell_idx_struct.tex2;

plot_trial_avg_imagesc(cell_struct,trial_indices,photo_ensembles,plot_cell_idx,plot_photo_types,plot_condition_types,opt)
export_fig([fig_save_path filesep 'TrialAvgStim_STATrace_' strrep(caiman_file,'.mat','.png')])

%%  == DIMENSIONALITY REDUCTION ==
%% params for factor analysis (also used for stim_opt so run this)
dr_opt.bin_size = 1;
dr_opt.gocue_bin = floor(opt.sta_gocue_frame/dr_opt.bin_size);
dr_opt.stim_bin = ceil(opt.sta_stim_frame/dr_opt.bin_size);
dr_opt.frame_rate = 30/dr_opt.bin_size;
dr_opt.Fs = opt.frame_rate;
dr_opt.trial_length = opt.trial_length/dr_opt.bin_size;
dr_opt.trial_color = trial_color;
dr_opt.idx_fields = {opt.trigger_idx_fd};

% fa_opt.avg_frames = fa_opt.stim_bin+30:1:fa_opt.gocue_bin;
dr_opt.avg_frames  = [];
easy_trial_idx = 1:10;

% dpca
% shuffle takes time - try parfor?
if strcmp(opt.method,'dpca')
dpca_opt = opt;
dpca_opt.trial_color = trial_color;
dr_opt.fd_names ={'stim_1_correct','stim_1_incorrect','stim_2_correct','stim_2_incorrect'};

[dpca_struct,traces_idx_struct,fa_trial_idx] = get_dpca_traj_brief(cell_struct,cell_idx_struct.(opt.trigger_idx_fd),dr_opt.fd_names,opt);

plot_pop_vectors(dpca_struct.traj_struct,dr_opt.fd_names,5,dpca_opt,...
    'plot_ylabel','PC level','plot_num_cols',2);
end
% run factor analysis
if strcmp(opt.method,'fa')
    % if FLAG_GOT_CHOICE_AUC
    dr_opt.fd_names ={'stim_1_correct','stim_1_incorrect','stim_2_correct','stim_2_incorrect'};
    % else
    %     fa_opt.fd_names ={'stim_1_correct','stim_2_correct'};
    %     disp('using correct trials only')
    % end
    dr_opt.plot_fds = dr_opt.fd_names;
    dr_opt.m = 3;
    dr_opt.IF_MEDFILT = 0;
    disp('Rnning factor analysis...')
    tic
    % excluding first 10 trials that were used for baselining - didnt make a
    % difference..
    [~,trial_idx_struct ]= cellfun(@(x)setdiff(trial_indices.(x),easy_trial_idx),dr_opt.fd_names,'UniformOutput',false);
    
    [traces_in,fa_trial_idx,num_trials] = get_input_seq(cell_struct,cell_idx_struct.(dr_opt.idx_fields{1}),...
        dr_opt.fd_names,dr_opt.bin_size,'IF_MEDFILT',dr_opt.IF_MEDFILT,...
        'fd_trial_idx',trial_idx_struct,'avg_frames',dr_opt.avg_frames);%%
    fa_struct = struct();
    fa_struct.mean = mean(traces_in);
    fa_struct.std = std(traces_in);
    [fa_struct.lambda,fa_struct.psi,fa_struct.T,fa_struct.stats,fa_struct.F] = factoran(traces_in,dr_opt.m,'Xtype','data','Maxit',1000);
    invsqrtPsi = diag(1 ./  sqrt(fa_struct.psi)); % get transition matrix (multiply to zscored data)
    fa_struct.transmat = invsqrtPsi/(fa_struct.lambda'*invsqrtPsi);
    if isempty(dr_opt.avg_frames)
        [fa_traj_struct,num_fs] = get_pop_vectors(fa_struct.F,dr_opt.trial_length,fa_trial_idx); % ran FA on full traces
    else
        full_traces_in  = get_input_seq(cell_struct,cell_idx_struct.(dr_opt.idx_fields{1}),...
            dr_opt.fd_names,dr_opt.bin_size,'IF_MEDFILT',dr_opt.IF_MEDFILT,...
            'fd_trial_idx',trial_idx_struct);%%
        F = (full_traces_in-fa_struct.mean)./fa_struct.std;
        [fa_traj_struct,num_fs] = get_pop_vectors(F,dr_opt.trial_length,fa_trial_idx);
    end
    plot_pop_vectors(fa_traj_struct,dr_opt.plot_fds,dr_opt.m,dr_opt,...
        'plot_ylabel','Factor level','plot_num_cols',2);
    toc
    disp('...Done')
end
    
% coding direction
% using raw onine traces to mimic online condition
if strcmp(opt.method,'cd')
choice_fds = {'stim_1_correct','stim_2_correct'};
cd_opt = opt;
cd_opt.trial_color = trial_color;
coding_direct = cell2mat(arrayfun(@(x)nanmean(mean(cell_struct(x).(choice_fds{1})(opt.sta_avg_frames,:),2)),cell_idx_struct.(opt.trigger_idx_fd),'UniformOutput', false))...
    - cell2mat(arrayfun(@(x)nanmean(mean(cell_struct(x).(choice_fds{2})(opt.sta_avg_frames,:),2)),cell_idx_struct.(opt.trigger_idx_fd),'UniformOutput', false));
% project to coding direction
cd_traj_struct = get_cd_projection(cell_struct(cell_idx_struct.(opt.trigger_idx_fd)),coding_direct,choice_fds);
plot_pop_vectors(cd_traj_struct,choice_fds,1,cd_opt,...
        'plot_ylabel','Projection')
end
%% == GET  DECODERS ==    
%% choose which trajectory to use 
switch opt.method 
    case 'fa'
    traj_struct = fa_traj_struct;
    case 'dpca'
    traj_struct = dpca_struct.traj_struct;
    case 'cd'
    traj_struct = cd_traj_struct;
end
% get stim decoder 
stim_opt = dr_opt;
stim_opt.fd_names = {'stim_1_correct','stim_2_correct'}; % stim 1 will be positive
stim_opt.frames_to_avg = stim_opt.stim_bin+30:1:stim_opt.gocue_bin;
stim_opt.frames_to_train = round([1:1:150]/stim_opt.bin_size);
stim_opt.Nstd = 2;
stim_opt.min_frames = 15;
stim_opt.IF_FRAMEWISE =0;
stim_struct = {};
stim_proj_struct = {};
disp('Running stim decoder...')
tic
stim_struct =  get_binary_classifier( stim_struct,traj_struct, stim_opt,...
    'IF_CROSSVAL',1,'IF_FRAMEWISE',stim_opt.IF_FRAMEWISE,'fd_names',stim_opt.fd_names);
 [stim_proj_struct] = get_projections(traj_struct,stim_struct.B(:,2:end)',stim_opt.fd_names,'proj_struct',stim_proj_struct,'bias',stim_struct.B(:,1));
plot_pop_vectors(stim_proj_struct,stim_opt.fd_names,1,stim_opt,...
        'plot_ylabel','Stim projection')
[ stim_struct ] =  get_binary_decoder_disc_time( stim_proj_struct, stim_struct,...
    stim_opt.fd_names,stim_opt,'IF_FRAMEWISE',stim_opt.IF_FRAMEWISE,'threshold',0);
toc
disp('Done')

figure; 
plot_binary_decoder(stim_struct,stim_opt)
suptitle('Stimulus decoder')
% export_fig([fig_save_path filesep 'StimDecoderPerform_' opt.save_name_ext strrep(online_file,'.mat','.png')])

%% get choice decoder (directly from traj, using threshold trials) 
choice_opt = dr_opt;
traj_struct.('go_trials') = cat(1,traj_struct.stim_1_incorrect,traj_struct.stim_2_correct);
traj_struct.('nogo_trials') = cat(1,traj_struct.stim_1_correct,traj_struct.stim_2_incorrect);

choice_opt.fd_names = {'stim_1_correct','stim_1_incorrect'}; %correct choice will be positive
choice_opt.fd_names = {'nogo_trials','go_trials'};
choice_opt.frames_to_avg = [120:150];
choice_opt.frames_to_train = round([1:1:150]/choice_opt.bin_size);
choice_opt.min_frames = 10;
choice_opt.Nstd = 1.5;
choice_opt.IF_FRAMEWISE = 0;
choice_struct = {};
choice_proj_struct = {};
disp('Rnning choice decoder...')
tic
choice_struct =  get_binary_classifier( choice_struct,traj_struct, choice_opt,...
    'IF_CROSSVAL',1,'IF_FRAMEWISE',choice_opt.IF_FRAMEWISE,'fd_names',choice_opt.fd_names);
 [choice_proj_struct] = get_projections(traj_struct,choice_struct.B(:,2:end)',choice_opt.fd_names,'proj_struct',choice_proj_struct,'bias',choice_struct.B(:,1));
plot_pop_vectors(choice_proj_struct,choice_opt.fd_names,1,choice_opt,...
        'plot_ylabel','Choice projection')
[ choice_struct ] =  get_binary_decoder_disc_time( choice_proj_struct, choice_struct,...
    choice_opt.fd_names,choice_opt,'IF_FRAMEWISE',choice_opt.IF_FRAMEWISE,'threshold',0);
toc
disp('Done')

figure; 
hold on
plot_binary_decoder(choice_struct,choice_opt)
suptitle('Choice decoder')
% export_fig([fig_save_path filesep 'ChoiceDecoderPerform_' opt.save_name_ext strrep(online_file,'.mat','.png')])

%% HMM-Mar
hmm_opt = dr_opt;
hmm_opt.plot_fds = hmm_opt.fd_names;
hmm_opt.num_states = 4;
hmm_opt.Fs = opt.frame_rate/hmm_opt.bin_size;            % sampling frequency
hmm_opt.K = hmm_opt.num_states;  	         % The number of states to infer
hmm_opt.order = 0; 	         % The lag used, this is only relevant when using MAR observations
hmm_opt.zeromean = 0; 	     % We do not want to model the mean, so zeromean is set on
hmm_opt.covtype = 'uniquefull';    % We want to model the full covariance matrix
hmm_opt.embeddedlags = 0; % 15 lags are used from -7 to 7
hmm_opt.pca = 0;   % The PCA dimensionality reduction is 2 times the number of ROIs because the covar matrix bins x channels is high
hmm_opt.initrep = 1; % to make it quicker - leave by default otherwise
hmm_opt.initcyc = 1; % to make it quicker - leave by default otherwise
hmm_opt.cyc = 200; % to make it quicker - leave by default otherwise
hmm_opt.verbose = 1;
hmm_opt.standardise = 1;
hmm_opt.max_num_run = 10; % maximum number of training cycle before states are related with trial structure
hmm_opt.avg_frames = [-10:1:0]+hmm_opt.gocue_bin; % use these frames to define state of a trial
hmm_opt.IF_MEDFILT = 0; % if filter input trace
hmm_opt.cell_idx_fd = 'tex';
hmm_opt.use_dpca = true;
disp('Rnning HMM analysis...')
tic
num_trials = round(size(dpca_struct.F,1)/hmm_opt.trial_length);
T = hmm_opt.trial_length * ones(num_trials,1);
if hmm_opt.use_dpca
    hmm_input = dpca_struct.F;
else
    hmm_input = get_input_seq(cell_struct,cell_idx_struct.(hmm_opt.cell_idx_fd),...
        hmm_opt.fd_names,hmm_opt.bin_size,'IF_MEDFILT',hmm_opt.IF_MEDFILT,'trial_dim',2,'time_dim',1);
end
[hmm,Gamma_hmm, ~, vpath] = hmmmar(hmm_input,T,hmm_opt);
toc
disp('...Done')
[hmm_struct,num_states] = get_pop_vectors(Gamma_hmm,hmm_opt.trial_length - hmm_opt.order,fa_trial_idx);
hmm_opt.m = num_states;
%plot states and paths
plot_pop_vectors(hmm_struct,hmm_opt.plot_fds,num_states,hmm_opt,'IF_MEDIAN',1,'plot_num_cols',4, 'plot_area',true);
vpath_struct = get_pop_vectors(vpath,hmm_opt.trial_length-hmm_opt.order,fa_trial_idx);
plot_pop_vectors(vpath_struct,hmm_opt.plot_fds,1,hmm_opt,...
    'ylimit',[0, hmm_opt.K+1],'plot_ylabel','Viterbi path','IF_MEDIAN',1,'plot_num_cols',4);
%% compare states during different trial epochs
figure('name','hmm state','position',[100 100 1200 400])
fd_colors =  cell2mat(cellfun(@(f)getfield(hmm_opt.trial_color,f)',hmm_opt.fd_names,'UniformOutput',false))';
subplot(1,4,1);
values = structfun(@(x)mode(x(:,1:opt.sta_baseline_frames),2),vpath_struct,'UniformOutput',false);
scatter_cmp_conditions(values,[],...
    1,fd_colors,'connect_scatter',0,'BriefXlabel',0,'ShowMeanInXlabel',1,'add_jitter',1,'plot_stats',1,'test_type','ttest');
ylabel('State index (baseline)')
xtickangle(45)

subplot(1,4,2);
values = structfun(@(x)mode(x(:,[-30:1:-1]+hmm_opt.gocue_bin),2),vpath_struct,'UniformOutput',false);
scatter_cmp_conditions(values,[],...
    1,fd_colors,'connect_scatter',0,'BriefXlabel',0,'ShowMeanInXlabel',1,'add_jitter',1,'plot_stats',1,'test_type','ttest');
ylabel('State index (withold)')
xtickangle(45)

subplot(1,4,3);
values = structfun(@(x)x(:,hmm_opt.gocue_bin),vpath_struct,'UniformOutput',false);
scatter_cmp_conditions(values,[],...
    1,fd_colors,'connect_scatter',0,'BriefXlabel',0,'ShowMeanInXlabel',1,'add_jitter',1,'plot_stats',1,'test_type','ttest');
ylabel('State index (go-cue)')
xtickangle(45)

subplot(1,4,4);
values = structfun(@(x)mode(x(:,[5:30]+hmm_opt.gocue_bin),2),vpath_struct,'UniformOutput',false);
scatter_cmp_conditions(values,[],...
    1,fd_colors,'connect_scatter',0,'BriefXlabel',0,'ShowMeanInXlabel',1,'add_jitter',1,'plot_stats',1,'test_type','ttest');
ylabel('State index (response)')
xtickangle(45)
%% probability and transition time to 'correct' state
% probability of being in the correct state at go-cue
probab = struct();
transtimes = struct();
onsets = getStateOnsets(vpath,T,hmm_opt.Fs,hmm.K);
[state_idx_struct,corr_struct,FLAG_MATCHED_STATES] = get_state_idx(hmm_struct,hmm_opt.avg_frames,hmm_opt.fd_names,num_states);

for f = 1:numel(hmm_opt.fd_names)
    fd = hmm_opt.fd_names{f};
    ref_fd = strrep(fd,'_incorrect','_correct');
    fdd = strrep(strrep(fd,'smooth_deconv_',''),'st_','');
    probab.(fdd) = hmm_struct.(fd)(:,hmm_opt.gocue_bin,state_idx_struct.(ref_fd));
    this_onsets = {onsets{fa_trial_idx.(fd),state_idx_struct.(ref_fd)}};
    this_onsets(cellfun(@(x)isempty(x),this_onsets))= {nan};
    this_onsets = cellfun(@(x)x(1),this_onsets);
    transtimes.(fdd) = this_onsets(this_onsets>1&this_onsets<opt.rw_win_end_sec);
end

figure('name','hmm summary','position',[100 100 1500 600])
this_title = [opt.exp_name, '_HmmProbTime'];
suptitle(strrep(this_title,'_',' '))
plot_color = cell2mat(cellfun(@(x)dr_opt.trial_color.(x)',hmm_opt.fd_names,'UniformOutput',false))';

subplot(1,2,1)
scatter_cmp_conditions(probab,[],...
    1,plot_color,'connect_scatter',0,'BriefXlabel',0,...
    'ShowMeanInXlabel',1,'add_jitter',1,'plot_stats',1);
ylabel('Probability (in correct state)')
subplot(1,2,2)
scatter_cmp_conditions(transtimes,[],...
    1,plot_color,'connect_scatter',0,'BriefXlabel',0,...
    'ShowMeanInXlabel',1,'add_jitter',1,'plot_stats',1);
ylabel('Transition time (into correct state)')

% export_fig([opt.save_path filesep this_title '.png'])

%% get decoder weights
stim_thresh = -stim_struct.B(:,1);
choice_thresh = -choice_struct.B(:,1);

switch opt.method
    case 'fa'
        stim_weights = fa_struct.transmat* stim_struct.B(:,2:end)';
        [stim_norm_weights,stim_norm_thresh] = get_norm_weights(stim_weights,stim_thresh,fa_struct.mean,fa_struct.std); % called in CompareDecoderWithAnimal
        choice_weights = fa_struct.transmat* choice_struct.B(:,2:end)';
        [choice_norm_weights,choice_norm_thresh] = get_norm_weights(choice_weights,choice_thresh,fa_struct.mean,fa_struct.std); % called in CompareDecoderWithAnimal
        
    case 'dpca'
        stim_weights = dpca_struct.W* stim_struct.B(:,2:end)';
        [stim_norm_weights,stim_norm_thresh] = get_norm_weights(stim_weights,stim_thresh,dpca_struct.mean,ones(size(dpca_struct.mean))); % called in CompareDecoderWithAnimal
        choice_weights = dpca_struct.W* choice_struct.B(:,2:end)';
        [choice_norm_weights,choice_norm_thresh] = get_norm_weights(choice_weights,choice_thresh,dpca_struct.mean,ones(size(dpca_struct.mean))); % called in CompareDecoderWithAnimal
        
        
    case 'cd'
        stim_weights =coding_direct'* stim_struct.B(2:end)';
        stim_norm_weights = stim_weights;
        stim_norm_thresh = stim_thresh;
        choice_norm_weights = coding_direct'* choice_struct.B(:,2:end)';
        choice_norm_thresh = choice_thresh;
        
end

stim_struct.norm_weights = stim_norm_weights;
stim_struct.norm_thresh = stim_norm_thresh;

choice_struct.norm_weights = choice_norm_weights;
choice_struct.norm_thresh = choice_norm_thresh;

%% Save decoders to structure
save_decod_struct = struct();
save_decod_struct.stim = stim_struct;
save_decod_struct.choice = choice_struct;
save_decod_struct.dpca_struct = dpca_struct;
output_save_name = [save_path filesep  'PosthocDecoders_' caiman_file ];
save(output_save_name,'save_decod_struct')
disp(['Output struct saved as:' output_save_name ])

%% compare decoder weights, color by stim auc
aucs = arrayfun(@(x)x.('correct_stimAUC_zscore'),cell_struct(cell_idx_struct.(opt.trigger_idx_fd)));
[auc_colors, sort_auc_colors ]= get_b2r_colors(aucs);
figure;scatter(choice_norm_weights,stim_norm_weights,[],auc_colors,'filled')
xlabel('Choice weights')
ylabel('Stim weights')
axis square
suptitle('Decoder weights')
colormap(sort_auc_colors)
h = colorbar;
ylabel(h,'Correct AUC')
set(h,'YTick',sort(aucs));
%% get projections on stim decoder
test_fd_names = {'stim_1_correct','stim_1_incorrect',...
                 'stim_2_correct','stim_2_incorrect',...
};

stim_proj_struct = struct();
[stim_proj_struct] = get_projections(cell_struct(cell_idx_struct.(opt.trigger_idx_fd)),stim_norm_weights,test_fd_names,'proj_struct',stim_proj_struct,'bias',-stim_norm_thresh,'IS_CELL_STRUCT',1);

plot_pop_vectors(stim_proj_struct,test_fd_names,1,opt,...
        'plot_ylabel','Projection','plot_num_cols',2,'IF_PLOT_RAW_ONLY',0)
suptitle('Stim decoder projections')

% compare projection on stim axis trials
figure('name','stim decoder projection','position',[100 100 800 400])
values = structfun(@(x)mean(x(:,opt.sta_avg_frames),2),stim_proj_struct,'UniformOutput',false);
fd_colors =  cell2mat(cellfun(@(f)getfield(trial_color,f),fields(values),'UniformOutput',false));
scatter_cmp_conditions(values,[],...
    1,fd_colors,'connect_scatter',0,'BriefXlabel',0,'ShowMeanInXlabel',1,'add_jitter',1);
suptitle('Stim decoder projections')
% export_fig([fig_save_path filesep 'StimDecoderProj_' opt.save_name_ext strrep(online_file,'.mat','.png')])


%% get projections on choice decoder
test_fd_names = {'stim_1_correct','stim_1_incorrect',...
                 'stim_2_correct','stim_2_incorrect',...
};
choice_proj_struct = struct();
[choice_proj_struct] = get_projections(cell_struct(cell_idx_struct.(opt.trigger_idx_fd)),choice_norm_weights,...
    test_fd_names,'proj_struct',choice_proj_struct,'bias',-choice_norm_thresh,'IS_CELL_STRUCT',1);
% compare projection on choice axis trials
plot_pop_vectors(choice_proj_struct,test_fd_names,1,opt,...
    'plot_ylabel','Projection','plot_num_cols',2,'IF_PLOT_RAW_ONLY',0)
suptitle('Choice decoder projections')

figure('name','choice decoder projection','position',[100 100 800 400])
values = structfun(@(x)mean(x(:,opt.sta_avg_frames),2),choice_proj_struct,'UniformOutput',false);
fd_colors =  cell2mat(cellfun(@(f)getfield(trial_color,f),fields(values),'UniformOutput',false));
scatter_cmp_conditions(values,[],...
    1,fd_colors,'connect_scatter',0,'BriefXlabel',0,'ShowMeanInXlabel',1,'add_jitter',1);
suptitle('Choice decoder projections')
% export_fig([fig_save_path filesep 'ChoiceDecoderProj_' strrep(online_file,'.mat','.png')])



%% plot accuracy - select tiral types to condition according to these plots
% compare decoder weights, color by stim auc
aucs = arrayfun(@(x)x.('correct_stimAUC_zscore'),cell_struct(cell_idx_struct.(opt.trigger_idx_fd)));
[auc_colors, sort_auc_colors ]= get_b2r_colors(aucs);
figure;scatter(choice_norm_weights,stim_norm_weights,[],auc_colors,'filled')
xlabel('Choice weights')
ylabel('Stim weights')
axis square
suptitle('Decoder weights')
colormap(sort_auc_colors)
h = colorbar;
ylabel(h,'Correct AUC')
set(h,'YTick',sort(aucs));
export_fig([fig_save_path filesep 'DecoderWeights' strrep(caiman_file,'.mat','.png')])

% test decoder accuracy
test_opt = stim_opt;
test_opt.fd_names = {'stim_1_correct','stim_1_incorrect',...
    'stim_2_correct','stim_2_incorrect'};

test_decoder_peformance(stim_proj_struct,test_opt,'Stim')
test_decoder_peformance(choice_proj_struct,test_opt,'Choice')


%% project catch trials on decoder
test_fd_names = {'stim_5_miss','stim_5_lick'};
[catch_proj_struct] = get_projections(cell_struct(cell_idx_struct.(opt.trigger_idx_fd)),pop_weights,...
    test_fd_names,'bias',-pop_thresh,'IS_CELL_STRUCT',1);
% compare projection on choice axis trials
plot_pop_vectors(catch_proj_struct,test_fd_names,1,opt,...
    'plot_ylabel','Projection','plot_num_cols',2,'IF_PLOT_RAW_ONLY',1)
suptitle('Catch trials decoder projections')

%% final check cells of interest
figure('name','trigger targets check','position',[100 100 2400 800]);
ax = subplot(1,3,1);
plot_value_in_rois( cell_struct, 'correct_stimAUC_zscore',[256 256],ax,...
    'IF_NORM_PIX',0,'IF_CONTOUR',0,'show_cell_idx',cell_idx_struct.(opt.trigger_idx_fd));
set(gca,'Ydir','reverse')
title('Texture selectivity (auc zscore)')

% mark texture cells on FOV
mark_cells_on_fov(cell_struct,cell_idx_struct.(opt.trigger_idx_fd),[.5 .5 .5],'MarkerSize',300,'Linewidth',2)
unique_target_fds = unique(opt.target_idx_fd);
for f = 1:numel(unique_target_fds)
    fd = strsplit(unique_target_fds{f},'_');
    fdd = fd{end};
    mark_cells_on_fov(cell_struct,cell_idx_struct.(unique_target_fds{f}),trial_color.(fdd),'MarkerSize',500,'Linewidth',2)
end

ax = subplot(1,3,3);
plot_value_in_rois( cell_struct, 'cnn_prediction',[256 256],ax,...
    'IF_NORM_PIX',0,'IF_CONTOUR',0,'show_cell_idx',cell_idx_struct.(opt.trigger_idx_fd));
set(gca,'Ydir','reverse')
title('CNN prediction score')

% mark texture cells on FOV
mark_cells_on_fov(cell_struct,cell_idx_struct.(opt.trigger_idx_fd),[.5 .5 .5],'MarkerSize',300,'Linewidth',2)
unique_target_fds = unique(opt.target_idx_fd);
for f = 1:numel(unique_target_fds)
    fd = strsplit(unique_target_fds{f},'_');
    fdd = fd{end};
    mark_cells_on_fov(cell_struct,cell_idx_struct.(unique_target_fds{f}),trial_color.(fdd),'MarkerSize',500,'Linewidth',2)
end

ax = subplot(1,3,2);
plot_value_in_rois( cell_struct, 'photo_auc_zscore',[256 256],ax,...
    'IF_NORM_PIX',0,'IF_CONTOUR',0,'show_cell_idx',cell_idx_struct.(opt.trigger_idx_fd));
set(gca,'Ydir','reverse')
title('Photo response (auc zscore)')

% mark texture cells on FOV
mark_cells_on_fov(cell_struct,cell_idx_struct.(opt.trigger_idx_fd),[.5 .5 .5],'MarkerSize',300,'Linewidth',2)
unique_target_fds = unique(opt.target_idx_fd);
for f = 1:numel(unique_target_fds)
    fd = strsplit(unique_target_fds{f},'_');
    fdd = fd{end};
    mark_cells_on_fov(cell_struct,cell_idx_struct.(unique_target_fds{f}),trial_color.(fdd),'MarkerSize',500,'Linewidth',2)
end

export_fig([fig_save_path filesep 'ValFOV_' strrep(caiman_file,'.mat','.png')])
%% check trigger traces
%% check trigger and target cell traces
figure('name','online trigger cell traces'); hold on
plot_offset = 20;
cell_count = 1;

this_idx = cell_idx_struct.(opt.trigger_idx_fd);
for i = 1:length(this_idx)
    ii = this_idx(i);
    cell_count = cell_count+1;
    plot(cell_struct(ii).filtC+cell_count*plot_offset,'black');
    text(double(caiman_data.t_init),cell_count*plot_offset,['Cell ' num2str(ii) ', W ' num2str(pop_weights(i))], 'horizontalalignment','right', 'color','black')
    
end
set(gca,'ytick',[])
ylim([0 plot_offset*(cell_count+1)])
for i = 1:numel(sens_stim_frames)
    this_color = trial_color.(['stim' num2str(trials.stim_type(i) )]);
    plot([sens_stim_frames(i) sens_stim_frames(i)]+opt.stimon_frame,ylim,'color',this_color) % trial-on + roughly first touch
    plot([sens_stim_frames(i) sens_stim_frames(i)]+opt.gocue_frame,ylim,'color',this_color,'linestyle',':') % go-cue
end


xlim([caiman_data.t_init tot_frames])





