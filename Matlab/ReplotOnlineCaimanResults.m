% replot online caiman results

root_dir = 'D:\TextureData\data\gonogo\';
all_session_paths = {
    'cb288\20200318\',... % dpca trace has interp frames
    'cb288\20200319\',... % fov unhealthy, has many bright dots
    'cb291\20200313\',...
    'cb273\20200318\',...
    'cb291\20200312\',...
    'cb291\20200309\'...
    };
%% init color
[trial_color] = online_tex_init_color();

%% load online saved-out baseline session struct 
session_idx = 1;
session_path = [root_dir filesep all_session_paths{session_idx}];
fileList = dir(fullfile([session_path '\analysis_files'],'*.mat'));
fileList = fileList(cellfun(@(x)contains(x,'ProcTex_'),{fileList(:).name}));
caiman_path = [session_path filesep 'analysis_files'];
caiman_file = fileList(1).name;
proc_data = load(fullfile(caiman_path,caiman_file));
tex_output_struct = proc_data.tex_output_struct;
clear proc_data
disp(['Loaded file :',fullfile(caiman_path,caiman_file)])
cell_struct = tex_output_struct.cell_struct;
cell_idx_struct = tex_output_struct.cell_idx_struct;
opt = tex_output_struct.opt;
opt.go_stim_types = 2;
opt.nogo_stim_types = 1;
opt.rw_win_end_sec = 5;
opt.withold_win_start_sec = 3;
opt.sta_avg_frames = [-30:1:-5]+opt.sta_gocue_frame; 
opt.IF_GO_NOGO = 1;

fig_save_path = [session_path filesep 'figures'];
opt.save_path = fig_save_path;
opt.trial_color = trial_color;
opt.sta_avg_frames = [120:130]; % overide to be before first photostim
%% add corrdinate to cells struct
cnm_dims = [256 256];
num_cells = size(cell_struct,2);
for i = 1:num_cells
    temp_coords = cell_struct(i).jsf.coordinates;
    lin_idx = zeros(size(temp_coords,1),1);
   
    for t = 1:size(temp_coords,1)
        lin_idx(t) = sub2ind(cnm_dims,temp_coords(t,1),temp_coords(t,2));
    end
    cell_struct(i).lin_coords = lin_idx;
    cell_struct(i).coordinates = cell_struct(i).jsf.coordinates;
    cell_struct(i).pix_values = cell_struct(i).jsf.values;
    cell_struct(i).centroid = cell_struct(i).jsf.centroid;
    cell_struct(i).weight = 0;
end
trig_cell_idx = cell_idx_struct.(opt.trigger_idx_fd);
trig_weights =  tex_output_struct.pop_params.weights;
for i = 1:numel(trig_cell_idx)
    t = trig_cell_idx(i);
    cell_struct(t).weight =trig_weights(i);
end
%% load baseline session caiman result
[caiman_file,caiman_path] = uigetfile('*.mat','Select baseline session caiman data',session_path);
caiman_data = load(fullfile(caiman_path,caiman_file));
disp(['Loaded file :',fullfile(caiman_path,caiman_file)])
cnm_accepted_idx = caiman_data.accepted_idx+1;
[cnm_struct,caiman_cell_struct] = get_cell_struct(caiman_data,[],[]);

%% load control session caiman result
[ctr_caiman_file,caiman_path] = uigetfile('*.mat','Select baseline control caiman data',session_path);
ctr_caiman_data = load(fullfile(caiman_path,ctr_caiman_file));
disp(['Loaded file :',fullfile(caiman_path,ctr_caiman_file)])
ctr_sens_stim_frames = ctr_caiman_data.sensory_stim_frames+ctr_caiman_data.t_init+opt.gocue_frame;

%% get online traj from baseline
%% get trial indices (baseline)
% trial types for go-nogo condition
% pyrtaoi stim: [1 2 3 4 3 4 5 5 1 2]; % trialOrder in caiman_data
% pyrtaoi var : [1 1 1 1 1 1 1 1 2 2]; % var = 2 are dummy closed-loop trials
% pybehav stim: [1 2 3 3 4 4 3 4 1 2];
% texture:      [1 2 3 3 3 3 3 3 1 2]; % by pybehav
% reward:       [0 2 0 0 2 2 2 0 0 2]; % by pybehav, 2 is go, 0 is no-go
% target:       [1 2 1 2 1 2 0 0 1 2]; % by pyrtaoi
disp('Select baseline pybehav')
[baseline_trials,tot_num_trials] = replot_get_trials_from_paq(opt,caiman_path);
[baseline_trials,baseline_trial_indices,sens_stim_frames] = sort_trial_types_baseline(baseline_trials,caiman_data,1,1,opt);

%exclude odd trials if any
disp('sorted trial indices')

%% load whisking suite2p animal struct to get whisking trials
animal_struct= load([session_path filesep 'analysis_files\baseline_animal_struct_engage']);
animal_struct = animal_struct.animal_struct;
baseline_whisk_trial_indices = animal_struct.whisk_trial_indices;

%% get trial indices (control session)
disp('Select control pybehav')
[ctr_trials,ctr_tot_num_trials] = replot_get_trials_from_paq(opt,caiman_path);
ctr_trials.trialOrder = ctr_caiman_data.trialOrder(1:ctr_tot_num_trials);
[ctr_photostim_trial_idx,caiman_num_photo_per_trial] = get_trials_with_photostim( ctr_caiman_data.sensory_stim_frames, ctr_caiman_data.online_photo_frames );
ctr_photostim_trial_idx(ctr_photostim_trial_idx>ctr_tot_num_trials) = [];
ctr_trials.photostim = zeros(size(ctr_trials.stim_type));
ctr_trials.photostim(ctr_photostim_trial_idx) = 1;
ctr_trials.trialVar=ones(1,ctr_tot_num_trials);
[ctr_trial_indices] = sort_trial_types_control(ctr_trials,opt);
%% load condition session caiman result
[caiman_file,caiman_path] = uigetfile('*.mat','Select codition session caiman data',caiman_path);
condition_caiman_data = load(fullfile(caiman_path,caiman_file));
disp(['Loaded file :',fullfile(caiman_path,caiman_file)])
%% load pybehav data (condition session)
disp('Select condition pybehav')
[trials,tot_num_trials] = replot_get_trials_from_paq(opt,caiman_path);
trials.trialOrder = condition_caiman_data.trialOrder(1:tot_num_trials);
try
trials.trialVar= condition_caiman_data.trialVar(1:tot_num_trials);
catch
    trials.trialVar = ones(1,tot_num_trials);
end
[photostim_trial_idx,num_photo_per_trial,first_photo_frames] = get_trials_with_photostim(  condition_caiman_data.sensory_stim_frames, condition_caiman_data.online_photo_frames );
[photostim_trial_idx1,num_photo_per_trial1,first_photo_frames1] = get_trials_with_photostim(  condition_caiman_data.stim_frames_caiman, condition_caiman_data.photo_stim_frames_caiman );
trials.photostim = zeros(1,tot_num_trials);
trials.photostim(photostim_trial_idx) = 1;
trials.photostim =trials.photostim (1:tot_num_trials);
first_photo_frames = first_photo_frames(1:tot_num_trials);
first_photo_frames1 = first_photo_frames1(1:tot_num_trials);
photostim_trial_idx = photostim_trial_idx(photostim_trial_idx<tot_num_trials);
[trial_indices] = sort_trial_types_condition(trials,opt);

% get online trial indices (include early licks)
trial_indices.stim_1_photo_online =intersect( photostim_trial_idx,find(condition_caiman_data.trialOrder==1));
trial_indices.stim_2_photo_online =intersect( photostim_trial_idx,find(condition_caiman_data.trialOrder==2));
trial_indices.stim_1_nonphoto_online =setdiff(find(condition_caiman_data.trialOrder==1),[photostim_trial_idx,1:10]);
trial_indices.stim_2_nonphoto_online =setdiff(find(condition_caiman_data.trialOrder==2), [photostim_trial_idx,1:10]);
trial_indices = structfun(@(x)x(x<tot_num_trials),trial_indices,'un',0);
condition_sens_stim_frames = condition_caiman_data.sensory_stim_frames+condition_caiman_data.t_init+opt.gocue_frame;

%% get lick times for each session (matching whisking sessions)
session_lick_times = {};
session_lick_times{1} = baseline_trials.firstlick;
session_lick_times{2} = trials.firstlick;
session_lick_times{3} = ctr_trials.firstlick;


%% load whisking kinematics and get kinematic sta
% 'amplitude' is the 'whisker envelope' used in Chen 2013
hs_path =  [session_path, filesep 'HsCam'];
hs_file = 'whisking_kinematics.mat';
whisking_data = load(fullfile(hs_path, hs_file)); % whisking
whisking = whisking_data.whisking;
clear whisking_data;
disp(['Loaded file :',fullfile(hs_path,hs_file)])
whisk_sta_struct = struct();
whisk_raw_sta_struct = struct();
fd = 'amplitude';
opt.whisk_trial_length = 1200;
opt.whisk_sta_gocue_frame = 400;
opt.whisk_frame_rate = 200;
num_sessions = numel(whisking);
whisk_amp_structs = {};
for s = 1:num_sessions
    
    this_raw_trace = whisking{s}.(fd);
    if strcmp(all_session_paths{session_idx},'cb291\202003007\')
        this_raw_trace = whisking{2}.(fd);
    end
    this_num_trials = floor(length(this_raw_trace)/opt.whisk_trial_length);
    this_sta_traces = reshape( this_raw_trace(1:this_num_trials*opt.whisk_trial_length),[opt.whisk_trial_length,this_num_trials])';
    whisk_sta_struct(s).raw_sta_traces = this_sta_traces;
    whisk_sta_struct(s).name = fd;
    % cound number of non-whisking bins in the trial
    all_lick_times = session_lick_times{s};
    this_num_trials = length(all_lick_times);
    all_whisk_amps = this_raw_trace;
    whisk_binsize = 20;
    whisk_amp_struct = this_sta_traces(1:this_num_trials,:);
    
    [all_trials_whisk_amp] = get_whisking_amount(whisk_amp_struct,all_lick_times,all_whisk_amps,whisk_binsize,opt);
    all_trials_whisk_amp = 1-all_trials_whisk_amp; 
    session_whisk_amp{s} = all_trials_whisk_amp;
    whisk_amp_structs{s} = whisk_amp_struct;
end

whisk_thresh = 0.5;
baseline_whisk_trial_indices = structfun(@(x)setdiff(x,find(session_whisk_amp{1}<whisk_thresh)),baseline_trial_indices,'un',0);
whisk_trial_indices = structfun(@(x)setdiff(x,find(session_whisk_amp{2}<whisk_thresh)),trial_indices,'un',0);
ctr_whisk_trial_indices = structfun(@(x)setdiff(x,find(session_whisk_amp{3}<whisk_thresh)),ctr_trial_indices,'un',0);
%% check lick time
trial_types = {'stim_1_correct','stim_1_incorrect','stim_2_correct','stim_2_incorrect'};
baseline_lick_struct = struct();
for t = trial_types
    baseline_lick_struct.(t{:}) = session_lick_times{1}(baseline_whisk_trial_indices.(t{:}));
end

%% get cell STA (baseline and control)
baseline_cell_struct = struct();
ctr_cell_struct = struct();
condition_cell_struct = struct();
% control 
skip_frames = ctr_caiman_data.frames_skipped + ctr_caiman_data.t_init;
ctr_num_frames = ctr_caiman_data.t_cnm;
ctr_tot_frames = ctr_num_frames + numel(skip_frames);
ctr_caiman_frames = setdiff([1:ctr_tot_frames],skip_frames);

% condition
num_frames = condition_caiman_data.t_cnm;
skip_frames = condition_caiman_data.frames_skipped + condition_caiman_data.t_init+1;
tot_frames = num_frames + numel(skip_frames);
caiman_frames = setdiff([1:tot_frames],skip_frames);


for i = 1:num_cells
    ii = cell_struct(i).cnm_idx;
    baseline_cell_struct(i).filtC =  caiman_data.filt_C(i,1:caiman_data.t_cnm); % filt_C only saved accepted cells
    baseline_cell_struct(i).deconvC =  caiman_data.cnm_C(ii,1:caiman_data.t_cnm); % filt_C only saved accepted cells
    baseline_cell_struct(i).cnm_idx = ii;
    
    temp_trace = nan(1,ctr_tot_frames);
    temp_trace(ctr_caiman_frames) =  ctr_caiman_data.filt_C(i, 1:ctr_num_frames);
    temp_trace = fillmissing(temp_trace,'linear');
    ctr_cell_struct(i).filtC =  temp_trace;
    
    temp_trace = nan(1,tot_frames);
    temp_trace(caiman_frames) =  condition_caiman_data.filt_C(i, 1:num_frames);
    temp_trace = fillmissing(temp_trace,'linear');
    condition_cell_struct(i).filtC =  temp_trace;
    
    temp_trace = nan(1,tot_frames);
    temp_trace(caiman_frames) =  condition_caiman_data.filt_C(i, 1:num_frames);
    temp_trace = fillmissing(temp_trace,'linear');
    condition_cell_struct(i).filtC =  temp_trace;

    
end


this_sens_stim_frames =   sens_stim_frames+opt.gocue_frame;
[baseline_cell_struct] = replot_get_cell_sta(baseline_cell_struct,this_sens_stim_frames,opt);
[ctr_cell_struct] = replot_get_cell_sta(ctr_cell_struct,ctr_sens_stim_frames,opt);
[condition_cell_struct] = replot_get_cell_sta(condition_cell_struct,condition_sens_stim_frames,opt);

disp('got cell_struct sta')

%% Normlalise traces to baseline (from first trials)
% get baseline and std from the first sets of easy trials
% - looks more consistent with traninig session to normalise this way
% doesnt make sense when sd is close to zero!! - just normalise to baseline
% in case intensity drifts...
easy_trial_idx = 1:10;
cell_bs = cell2mat(arrayfun(@(x)mean(reshape(baseline_cell_struct(x).raw_sta_traces(easy_trial_idx,1:opt.sta_baseline_frames),[1,numel(easy_trial_idx)*opt.sta_baseline_frames]),2),1:num_cells,'UniformOutput', false));
cell_sd = cell2mat(arrayfun(@(x)std(reshape(baseline_cell_struct(x).raw_sta_traces(easy_trial_idx,:),[1,numel(easy_trial_idx)*opt.trial_length])),1:num_cells,'UniformOutput', false));
cell_mean = cell2mat(arrayfun(@(x)mean(baseline_cell_struct(x).raw_sta_trace(:)),1:num_cells,'UniformOutput', false));

for i = 1:num_cells
    baseline_cell_struct(i).raw_sta_traces = (baseline_cell_struct(i).raw_sta_traces - cell_bs(i));
    baseline_cell_struct(i).raw_sta_trace = (baseline_cell_struct(i).raw_sta_trace - cell_bs(i));
    ctr_cell_struct(i).raw_sta_traces = ctr_cell_struct(i).raw_sta_traces-ctr_caiman_data.bs_level(i);
    condition_cell_struct(i).raw_sta_traces = condition_cell_struct(i).raw_sta_traces-condition_caiman_data.bs_level(i);
end

disp('normalised cell_struct raw_sta_traces')

%% Sort STAs by different trial types (baseline and control)
trial_types = fields(baseline_trial_indices);
baseline_raw_cell_struct = struct();

for i = 1:numel(trial_types)
    this_fd = trial_types{i};
    this_idx = baseline_trial_indices.(this_fd);
    this_idx = this_idx(this_idx<=length(baseline_trials.stim_type));
%     this_idx = setdiff(this_idx,easy_trial_idx);
    if ~isempty(this_idx)
        for c = 1:num_cells
            baseline_raw_cell_struct(c).(this_fd) = baseline_cell_struct(c).raw_sta_traces( this_idx,:)';
        end
    end
end

trial_types = fields(ctr_trial_indices);
num_sta_traces = size(ctr_cell_struct(1).raw_sta_traces,1);
ctr_raw_cell_struct = struct();

for i = 1:numel(trial_types)
    this_fd = trial_types{i};
    this_idx = ctr_trial_indices.(this_fd);
    this_idx = this_idx(this_idx<=min([ctr_tot_num_trials,num_sta_traces]));
%     this_idx = setdiff(this_idx,easy_trial_idx);
    if ~isempty(this_idx)
        for c = 1:num_cells
            ctr_raw_cell_struct(c).(this_fd) = ctr_cell_struct(c).raw_sta_traces( this_idx,:)';          
        end
    end
end

disp('sorted raw_sta_traces')
%% see if non-whisking trials can be exluded using background neuropil 
opt.bg_sta_avg_frames = 90:150;
background_sta_traces = struct();
background_whisk_sta_traces = struct();
whisk_amp_sta = struct(); 
backgroundC = zscore(caiman_data.noisyC(1,1:num_frames));
[~,~,~,bg_sta_traces] = make_sta_from_traces(backgroundC,this_sens_stim_frames,opt.sta_pre_frames,opt.sta_post_frames,1:opt.sta_baseline_frames);
% figure;plot(backgroundC);
trial_types = {'stim_1_correct','stim_1_incorrect','stim_2_correct','stim_2_incorrect'};
whisk_trial_types = trial_types;

for i = 1:numel(trial_types)
    this_fd = trial_types{i};
    this_idx = baseline_trial_indices.(this_fd);
    this_idx = this_idx(this_idx<=min([tot_num_trials]));
    background_sta_traces.(strrep(this_fd,'_nonphoto','')) = bg_sta_traces( this_idx,:)';
    
    whisk_amp_sta.(this_fd) = animal_struct.all_trials_whisk_amp(this_idx);
    

    this_fd = whisk_trial_types{i};
    this_idx = baseline_whisk_trial_indices.(this_fd);
    this_idx = this_idx(this_idx<=min([tot_num_trials]));
    background_whisk_sta_traces.(strrep(this_fd,'_nonphoto','')) = bg_sta_traces( this_idx,:)';
    
end

bg_trial_avg = structfun(@(x)mean(x(opt.bg_sta_avg_frames ,:),1),background_sta_traces,'un',0);

bg_whisk_trial_avg = structfun(@(x)mean(x(opt.bg_sta_avg_frames ,:),1),background_whisk_sta_traces,'un',0);

all_trial_bg = cell2mat(struct2cell(bg_trial_avg)');
all_trial_whisk_amp = cell2mat(struct2cell(whisk_amp_sta)');
whisk_trial_bg = cell2mat(struct2cell(bg_whisk_trial_avg)');

% bimodality test
[xpdf] = compute_xpdf(all_trial_bg);
[dip, p_dip] = HartigansDipSignifTest(xpdf, 1000);
[b] = bmtest(all_trial_bg);

[counts] = hist(all_trial_bg./(max(all_trial_bg)-min(all_trial_bg)),128);
level = otsuthresh(counts).*(max(all_trial_bg)-min(all_trial_bg))+min(all_trial_bg);


figure('name','neuropil vs whisking','position',[100 100 1200 400])
subplot(1,3,1)
binwidth = 0.02;
hold on
histogram(all_trial_bg,'displaystyle','stairs','edgecolor',[.5 .5 .5],'displayname','all trials','binwidth',binwidth,'linewidth',2)
histogram(whisk_trial_bg,'displaystyle','stairs','edgecolor','black','displayname','whisking trials','binwidth',binwidth,'linewidth',2)
plot(level.*[1,1],ylim,'color','r','linestyle',':','displayname','threshold');
box off
xlabel('Background neuropil avg.')
ylabel('Num. trials')
legend('location','northoutside','box','off')

subplot(1,3,2)
hold on
scatter(all_trial_bg,all_trial_whisk_amp,'markeredgecolor','black','markerfacecolor','none');
plot(level.*[1,1],ylim,'color','r','linestyle',':','displayname','threshold');
axis square
ylabel('whisking amp.')
xlabel('background calcium (neuropil)')

subplot(1,3,3)
fd_colors =  cell2mat(cellfun(@(f)getfield(trial_color,f),fields(bg_trial_avg),'UniformOutput',false));
scatter_cmp_conditions(bg_trial_avg,[],...
    1,fd_colors,'plot_stats',1,'connect_scatter',0,'BriefXlabel',0,'ShowMeanInXlabel',1,'add_jitter',1,'xtick_label',{'CR','FA','Hit','Miss'});
axis square
xtickangle(45)
ylabel('Background neuropil avg.')

% export_fig([fig_save_path filesep 'CnmBackgroundNeuropilVsWhiskingAmp' '.png'])
%% =========  Post-hoc Decoder using whiking trials only ===========
% 1. sta traces from trials with higher calcium (potentially high whisking)
baseline_train_trial_indices = struct();
for i = 1:numel(trial_types)
    this_fd = trial_types{i};  
    baseline_train_trial_indices.(this_fd) = baseline_trial_indices.(this_fd)(bg_trial_avg.(this_fd)>level);
end

trial_types = fields(baseline_train_trial_indices);
baseline_train_cell_struct = struct();% raw online traces (filtC)
baseline_train_deconv_struct = struct();% deconv traces
baseline_train_deconv_struct1 = struct();
baseline_train_deconv_struct2 = struct();
for i = 1:numel(trial_types)
    this_fd = trial_types{i};
    this_idx = baseline_train_trial_indices.(this_fd);
    this_idx = this_idx(this_idx<=min([tot_num_trials]));
%     this_idx = setdiff(this_idx,easy_trial_idx);
    if ~isempty(this_idx)
        for c = 1:num_cells
            baseline_train_cell_struct(c).(strrep(this_fd,'_nonphoto','')) = baseline_cell_struct(c).raw_sta_traces( this_idx,:)';
            baseline_train_deconv_struct(c).(strrep(this_fd,'_nonphoto','')) = baseline_cell_struct(c).sta_traces( this_idx,:)';
            
            % split into two halves
            this_num_trials = length(this_idx);
            baseline_train_deconv_struct1(c).(strrep(this_fd,'_nonphoto','')) = baseline_cell_struct(c).sta_traces( this_idx(1:ceil(this_num_trials/2)),:)';
            baseline_train_deconv_struct2(c).(strrep(this_fd,'_nonphoto','')) = baseline_cell_struct(c).sta_traces( this_idx(setdiff(1:this_num_trials,1:ceil(this_num_trials/2))),:)';

        end
    end
    baseline_train_trial_indices.(strrep(this_fd,'_nonphoto','')) = this_idx;
    baseline_train_bg_struct.(strrep(this_fd,'_nonphoto','')) = bg_sta_traces(this_idx,:);
end

%% 2. get cell identity
opt.N = 1.5;
opt.flag_use_peak = 0;
opt.sta_peak_search_range = [120:150];
[baseline_train_deconv_struct] = get_cell_auc(baseline_train_deconv_struct,{'stim_1_correct','stim_2_correct'},'correct_stimAUC_whisk',opt);
[baseline_train_deconv_struct] = get_cell_auc(baseline_train_deconv_struct,{'stim_1_correct','stim_1_incorrect'},'nogochoiceAUC_whisk',opt);
[baseline_train_deconv_struct] = get_cell_auc(baseline_train_deconv_struct,{'stim_2_incorrect','stim_2_correct'},'gochoiceAUC_whisk',opt);

[baseline_train_deconv_struct1] = get_cell_auc(baseline_train_deconv_struct1,{'stim_1_correct','stim_2_correct'},'correct_stimAUC_whisk',opt);
[baseline_train_deconv_struct2] = get_cell_auc(baseline_train_deconv_struct2,{'stim_1_correct','stim_2_correct'},'correct_stimAUC_whisk',opt);

[baseline_train_deconv_struct] = get_cell_auc(baseline_train_deconv_struct,{'stim_1_incorrect','stim_2_correct'},'stimAUC_fahit',opt);
[baseline_train_deconv_struct] = get_cell_auc(baseline_train_deconv_struct,{'stim_1_correct','stim_2_incorrect'},'stimAUC_crmiss',opt);
[baseline_train_cell_struct] = get_cell_active_auc(baseline_train_cell_struct,{'stim_1_correct','stim_2_correct'},opt);


correct_stimulusAUC_zscore = extractfield(baseline_train_deconv_struct,'correct_stimAUC_whisk_zscore');
gochoiceAUC_zscore = extractfield(baseline_train_deconv_struct,'gochoiceAUC_whisk_zscore');
nogochoiceAUC_zscore = extractfield(baseline_train_deconv_struct,'nogochoiceAUC_whisk_zscore');
correct_stimulusAUC_zscore1 = extractfield(baseline_train_deconv_struct1,'correct_stimAUC_whisk_zscore');
correct_stimulusAUC_zscore2 = extractfield(baseline_train_deconv_struct2,'correct_stimAUC_whisk_zscore');
lickstimAUC_zscore = extractfield(baseline_train_deconv_struct,'stimAUC_fahit_zscore');
nolickstimAUC_zscore = extractfield(baseline_train_deconv_struct,'stimAUC_crmiss_zscore');
active_stim1_auc = extractfield(baseline_train_cell_struct,'stim_1_correct_activeAUC');
active_stim2_auc = extractfield(baseline_train_cell_struct,'stim_2_correct_activeAUC');

cell_idx_struct.whisk_tex2 = find(correct_stimulusAUC_zscore>opt.N&active_stim2_auc>0.4); % cells prefering texture1 in correct trials
cell_idx_struct.whisk_tex1 = find(correct_stimulusAUC_zscore<-opt.N&active_stim1_auc>0.4); % cells prefering texture2 in correct trials
cell_idx_struct.whisk_tex = unique([cell_idx_struct.whisk_tex1, cell_idx_struct.whisk_tex2]);
cell_idx_struct.robust_tex2 = find(correct_stimulusAUC_zscore1>opt.N&correct_stimulusAUC_zscore2>opt.N&active_stim2_auc>0.4);
cell_idx_struct.robust_tex1 = find(correct_stimulusAUC_zscore1<-opt.N&correct_stimulusAUC_zscore2<-opt.N&active_stim1_auc>0.4);
cell_idx_struct.robust_tex = unique([cell_idx_struct.robust_tex1 cell_idx_struct.robust_tex2]);

% cells that robustly encode stim or choice
cell_idx_struct.go =  find(gochoiceAUC_zscore>opt.N & nogochoiceAUC_zscore>opt.N);
cell_idx_struct.nogo = find(gochoiceAUC_zscore<-opt.N & nogochoiceAUC_zscore<-opt.N);
cell_idx_struct.stim2 =  find(lickstimAUC_zscore>opt.N & nolickstimAUC_zscore>opt.N&active_stim2_auc>0.4);
cell_idx_struct.stim1 =  find(lickstimAUC_zscore<-opt.N & nolickstimAUC_zscore<-opt.N&active_stim1_auc>0.4);
cell_idx_struct.stim1_nogo = intersect(cell_idx_struct.stim1,cell_idx_struct.nogo);
cell_idx_struct.choice2 =intersect(cell_idx_struct.whisk_tex2, find(abs(gochoiceAUC_zscore)>opt.N&active_stim2_auc>0.4)); % exclude cells with higher activity in baseline than withold window
cell_idx_struct.choice1 =intersect(cell_idx_struct.whisk_tex1, find(abs(nogochoiceAUC_zscore)>opt.N&active_stim1_auc>0.4));

cell_idx_struct.choice2 =find(abs(gochoiceAUC_zscore)>opt.N&active_stim2_auc>0.4); % exclude cells with higher activity in baseline than withold window
cell_idx_struct.choice1 =find(abs(nogochoiceAUC_zscore)>opt.N&active_stim1_auc>0.4);

cell_idx_struct.choice = unique([cell_idx_struct.choice1   cell_idx_struct.choice2]);
cell_idx_struct.stim = unique([cell_idx_struct.stim1   cell_idx_struct.stim2]);
cell_idx_struct.stimchoice = unique([cell_idx_struct.stim,cell_idx_struct.choice]);

cell_idx_struct.stim = unique([cell_idx_struct.stim1 cell_idx_struct.stim2]);
cell_idx_struct.intersect = intersect(cell_idx_struct.stim, cell_idx_struct.choice);
cell_idx_struct.intersect_stim1 =  intersect(cell_idx_struct.stim1, cell_idx_struct.choice);
cell_idx_struct.nogo_stim2 =  intersect(cell_idx_struct.stim2, cell_idx_struct.nogo);
cell_idx_struct.go_stim2 =  intersect(cell_idx_struct.stim2, cell_idx_struct.go);
cell_idx_struct.nogo_stim1 =  intersect(cell_idx_struct.stim1, cell_idx_struct.nogo);
cell_idx_struct.go_stim1 =  intersect(cell_idx_struct.stim1, cell_idx_struct.go);
cell_idx_struct.go_only = setdiff( find(gochoiceAUC_zscore>opt.N),abs(lickstimAUC_zscore)<opt.N);

% stimulus encoding cells that does not have the opposite choice preference
cell_idx_struct.stim_exl = unique([setdiff(cell_idx_struct.stim1,cell_idx_struct.go) setdiff(cell_idx_struct.stim2,cell_idx_struct.nogo)]);
cell_idx_struct.stim2_exl = setdiff(cell_idx_struct.stim2,cell_idx_struct.nogo);


figure('name','choice auc')
mark_cell_types = {'stim1','stim2','choice1','choice2'};
mark_cell_colors = {'stim1','stim2','hold','lick',};
for dummy = 1
subplot(1,2,1)
hold on
xx = gochoiceAUC_zscore;
yy = nogochoiceAUC_zscore;
scatter(xx,yy,'markeredgecolor',[.7 .7 .7])

[ r,~,~,~,~,~,p ] = plot_fit_linear(xx,yy,[min(xx),max(xx)],[.6 .6 .6]);
text(1,1,['r=' num2str(r,'%10.2f') '; P= ' num2str(p,'%10.4f')],'units','normalized', 'horizontalalignment','right','verticalalignment','top', 'color','black')
xlim([-5 5])
ylim([-5 5])
plot([1,1].*opt.N,ylim,':','color','black')
plot([-1,-1].*opt.N,ylim,':','color','black')
plot(xlim,[-1,-1].*opt.N,':','color','black')
plot(xlim,[1,1].*opt.N,':','color','black')
plot([0,0].*opt.N,ylim,':','color','black')
plot(xlim,[0,0].*opt.N,':','color','black')
p = [];
for i = 1:numel(mark_cell_types)
    this_cell_idx = cell_idx_struct.(mark_cell_types{i});
    p(i)=scatter(xx(this_cell_idx),yy(this_cell_idx),'markerfacecolor',trial_color.(mark_cell_colors{i}),'markeredgecolor',trial_color.(mark_cell_colors{i}),...
        'markerfacealpha',.5,'markeredgealpha',.5);
end
axis square
xlabel('Stimulus selectivity (FA vs Hit)')
ylabel('Stimulus selectiviey (CR vs Miss)')
legend(p,mark_cell_types,'location','eastoutside','box','off')
axis square
xlabel('Choice selectivity (Miss vs Hit)')
ylabel('Choice selectiviey (CR vs FA)')


subplot(1,2,2)
hold on
xx = lickstimAUC_zscore;
yy = nolickstimAUC_zscore;
scatter(xx,yy,'markeredgecolor',[.6 .6 .6])

[ r,~,~,~,~,~,p ] = plot_fit_linear(xx,yy,[min(xx),max(xx)],[.6 .6 .6]);
text(1,1,['r=' num2str(r,'%10.2f') '; P= ' num2str(p,'%10.4f')],'units','normalized', 'horizontalalignment','right','verticalalignment','top', 'color','black')
xlim([-5 5])
ylim([-5 5])
plot([1,1].*opt.N,ylim,':','color','black')
plot([-1,-1].*opt.N,ylim,':','color','black')
plot(xlim,[-1,-1].*opt.N,':','color','black')
plot(xlim,[1,1].*opt.N,':','color','black')
plot([0,0].*opt.N,ylim,':','color','black')
plot(xlim,[0,0].*opt.N,':','color','black')
p = [];
for i = 1:numel(mark_cell_types)
    this_cell_idx = cell_idx_struct.(mark_cell_types{i});
    p(i)=scatter(xx(this_cell_idx),yy(this_cell_idx),'markerfacecolor',trial_color.(mark_cell_colors{i}),'markeredgecolor',[.6 .6 .6],'markerfacealpha',.5);
end
axis square
xlabel('Stimulus selectivity (FA vs Hit)')
ylabel('Stimulus selectiviey (CR vs Miss)')
legend(p,mark_cell_types,'location','eastoutside','box','off')
end

%% get ensemble sta traces
opt.sta_avg_frames = [120:150];
cell_types = fields(cell_idx_struct);
plot_fds = {'stim_1_correct','stim_1_incorrect','stim_2_correct','stim_2_incorrect'};
plot_fds_names = {'CR','FA','Hit','Miss'};
for f = 1:numel(cell_types)
    for ff = 1:numel(plot_fds)
        fd = plot_fds{ff};
        ensemble_mean.(cell_types{f}).(fd) = [];
        ensemble_traces.(cell_types{f}).(fd) = [];
    end
end

[ensemble_traces,ensemble_mean,ensemble_traces_raw,ensemble_mean_raw,ensemble_cell_mean] = get_ensemble_mean(ensemble_mean,ensemble_traces,baseline_train_deconv_struct,cell_idx_struct,plot_fds,opt);
whisk_avg_frames = round(((opt.sta_avg_frames - opt.sta_gocue_frame)/opt.frame_rate)*opt.whisk_frame_rate+opt.whisk_sta_gocue_frame);
baseline_whisk_amps = structfun(@(x)mean(whisk_amp_structs{1}(x,whisk_avg_frames),2),baseline_train_trial_indices,'un',0);


condition_whisk_amps = struct();
for p = 1:numel(plot_fds)
condition_whisk_amps.(plot_fds{p}) = mean(whisk_amp_structs{2}(whisk_trial_indices.(plot_fds{p}),whisk_avg_frames),2);
end
%% plot ensemble sta traces
figure('name','ensemble traces','units','normalized','outerposition',[0 0 1 .8])
plot_cell_types = {'all','stim1','stim2','choice1','choice2'};
num_plots = numel(plot_cell_types)+1;

trace_plot_fds = {'stim_1_correct','stim_1_incorrect','stim_2_correct','stim_2_incorrect'};

subplot(2,num_plots,1)
hold on
% whisking amp
this_xticks = ([1:opt.whisk_trial_length]-opt.whisk_sta_gocue_frame)./opt.whisk_frame_rate;

for f = 1:numel(trace_plot_fds)
    fd = trace_plot_fds{f};
    this_traces = whisk_amp_structs{1}(baseline_train_trial_indices.(fd),:);
    shadedErrorBar(this_xticks,mean(this_traces,1), std(this_traces,[],1),{'color',trial_color.(fd),'linewidth',2},0.1);
end
plot([0 0],ylim,':','color',[.5 .5 .5],'linewidth',2) % go-cue

ylabel([ 'Whisking amp.'])
xlabel('Time from go-cue (sec)')
xlim([-3,4])
axis square


subplot(2,num_plots,1+num_plots)
fd_colors =  cell2mat(cellfun(@(f)getfield(trial_color,f),fields(baseline_whisk_amps),'UniformOutput',false));
scatter_cmp_conditions(baseline_whisk_amps,[],...
    1,fd_colors,'connect_scatter',0,'BriefXlabel',1,'ShowMeanInXlabel',1,'add_jitter',1,...
    'plot_stats',1,'test_type','ranksum','xtick_label',plot_fds_names);
ylabel('Whisking amp.')
axis square

this_xticks = ([1:opt.trial_length]-opt.sta_gocue_frame)./opt.frame_rate;

plot_count = 2;
for p = 1:numel(plot_cell_types)
    subplot(2,num_plots,plot_count)
    hold on
    
    cell_id = plot_cell_types{p};
    this_num_cells = numel(cell_idx_struct.(cell_id));
    for f = 1:numel(trace_plot_fds)
        fd = trace_plot_fds{f};
        this_traces = ensemble_traces_raw.(cell_id).(fd)';
        shadedErrorBar(this_xticks,mean(this_traces,1), std(this_traces,[],1),{'color',trial_color.(fd),'linewidth',2},0.1);
    end
    title([ strrep(cell_id,'_','-') ' ens. avg., ' num2str(this_num_cells) ' cells' ])
    ylabel('Activity level')
    xlabel('Time from go-cue (sec)')
    ylim([-0.5,5])
    xlim([-5,4])

    plot([0 0],ylim,':','color',[.5 .5 .5],'linewidth',2) % go-cue
    plot([-1 0],[0 0],'color','black','linewidth',2) % sta average window (withold, 1sec before gocue)
    axis square
    
    subplot(2,num_plots,plot_count+num_plots)

    this_values = ensemble_mean_raw.(cell_id);
    fd_colors =  cell2mat(cellfun(@(f)getfield(trial_color,f),fields(this_values),'UniformOutput',false));
    
    scatter_cmp_conditions(this_values,[],...
        1,fd_colors,'connect_scatter',0,'BriefXlabel',1,'ShowMeanInXlabel',1,'add_jitter',1,...
        'plot_stats',1,'test_type','ranksum','xtick_label',plot_fds_names);
    ylabel([ strrep(cell_id,'_','-') ' cell avg.' ])
    axis square
     plot_count =plot_count+ 1;
   
end
% export_fig([fig_save_path filesep 'OnlineTexCellsEnsembleMean' '.png'])
%% plot time course of choice-related signal
figure('name','diff traces','units','normalized','outerposition',[0 0 1 1])
this_xticks = ([1:opt.trial_length]-opt.sta_gocue_frame)./opt.frame_rate;
plot_count = 1;
num_plots = numel(plot_cell_types);
withold_avg_frames = opt.sta_avg_frames;
baseline_avg_frames = 30:60; % when motor is rotating
withold_time = (withold_avg_frames-opt.sta_gocue_frame)./opt.frame_rate;
bs_time = (baseline_avg_frames-opt.sta_gocue_frame)./opt.frame_rate;

for p = 1:numel(plot_cell_types)
    cell_id = plot_cell_types{p};

    subplot(4,num_plots,plot_count)
    title([strrep(cell_id,'_',' ') 'cells'])
    hold on
    this_traces = ensemble_cell_mean.(cell_id).('stim_1_incorrect')' - ensemble_cell_mean.(cell_id).('stim_1_correct')';
    shadedErrorBar(this_xticks,mean(this_traces,1), std(this_traces,[],1),{'color',[.5 .5 .5],'linewidth',2},0.3);
%     plot(repmat(this_xticks,[size(this_traces,1),1])',this_traces','color',[.5 .5 .5])
%       shadedErrorBar(this_xticks,median(this_traces,1), abs(quantile(this_traces,[.25 .75])-median(this_traces,1)),{'color',[.5 .5 .5],'linewidth',2},0.1);
  
     xlim([-4,4])
%     ylim([-1,1])
    plot([0 0],ylim,':','color',[.5 .5 .5],'linewidth',2) % go-cue
    plot(xlim,[0 0],':','color','black','linewidth',1) 
    plot([bs_time(1) bs_time(end)],[0 0],'color',[.3 .3 .3],'linewidth',2) 

    plot([-1 0],[0 0],'color','black','linewidth',2) % sta average window (withold, 1sec before gocue)
    axis square
    ylabel([ 'Go-activity (FA-CR)' ])
    xlabel('Time from go-cue (sec)')

    subplot(4,num_plots,plot_count+num_plots)
    hold on
    this_diff = struct();
    this_diff.baseline = mean(this_traces(:,baseline_avg_frames),2);
    this_diff.withhold = mean(this_traces(:,withold_avg_frames),2);

    scatter_cmp_conditions(this_diff,[],...
        1,[0 0 0;0 0 0],'connect_scatter',1,'BriefXlabel',1,'ShowMeanInXlabel',0,'add_jitter',1,...
        'plot_stats',1,'test_type','signrank');
    xtickangle(45)
    ylabel(['Trial avg. FA-CR' ])
    xlabel([])
    xlim([0 3])
    axis square
    
    
    subplot(4,num_plots,plot_count+num_plots*2)
    title([strrep(cell_id,'_',' ') 'cells'])
    hold on
    this_traces = ensemble_cell_mean.(cell_id).('stim_2_correct')' - ensemble_cell_mean.(cell_id).('stim_2_incorrect')';
    shadedErrorBar(this_xticks,mean(this_traces,1), std(this_traces,[],1),{'color',[.5 .5 .5],'linewidth',2},0.3);
%     plot(repmat(this_xticks,[size(this_traces,1),1])',this_traces','color',[.5 .5 .5])
%       shadedErrorBar(this_xticks,median(this_traces,1), abs(quantile(this_traces,[.25 .75])-median(this_traces,1)),{'color',[.5 .5 .5],'linewidth',2},0.1);
  
     xlim([-4,4])
%     ylim([-1,1])
    plot([0 0],ylim,':','color',[.5 .5 .5],'linewidth',2) % go-cue
    plot(xlim,[0 0],':','color','black','linewidth',1) 
    plot([bs_time(1) bs_time(end)],[0 0],'color',[.3 .3 .3],'linewidth',2) 

    plot([-1 0],[0 0],'color','black','linewidth',2) % sta average window (withold, 1sec before gocue)
    axis square
    ylabel([ 'Go-activity (Hit-Miss)' ])
    xlabel('Time from go-cue (sec)')

    subplot(4,num_plots,plot_count+num_plots*3)
    hold on
    this_diff = struct();
    this_diff.baseline = mean(this_traces(:,baseline_avg_frames),2);
    this_diff.withhold = mean(this_traces(:,withold_avg_frames),2);

    scatter_cmp_conditions(this_diff,[],...
        1,[0 0 0;0 0 0],'connect_scatter',1,'BriefXlabel',1,'ShowMeanInXlabel',0,'add_jitter',1,...
        'plot_stats',1,'test_type','signrank');
    xtickangle(45)
    ylabel(['Trial avg. Hit-Miss' ])
    
    xlim([0 3])
     xlabel([])

    axis square

    plot_count =plot_count+ 1;
   
end

%% correlate ensemble activity with whisking, trialwise
figure('name','compare baseline','units','normalized','outerposition',[0 0 1 1])
for f = 1:numel(trace_plot_fds)
    fd = trace_plot_fds{f};
    whisk_amp_sta_struct.(fd) = whisk_amp_structs{1}(baseline_train_trial_indices.(fd),:);
end
[whisk_withold_amp,neural_withold_amp] = plot_cmp_neural_whisk_amp(whisk_amp_sta_struct,ensemble_traces_raw.('whisk_tex1'),whisk_avg_frames,trace_plot_fds,...
    opt.sta_avg_frames,'Whisking amp.','Ensemble mean',trial_color,plot_fds_names);


%% 3. get trajectory and choice decoder (directly from traj, using threshold trials) 
opt.method = 'dpca';
choice_opt = dpca_opt;
choice_opt.frames_to_avg = [130:145];
choice_opt.frames_to_train = round([1:1:150]/choice_opt.bin_size);
choice_opt.min_frames = 10;
choice_opt.Nstd = 1.5;
choice_opt.IF_FRAMEWISE = 0;
choice_opt.Fs = 30;

dpca_opt = opt;
dpca_opt.trial_color = trial_color;
dpca_opt.bin_size = 1;
dpca_opt.gocue_bin = 150;
dpca_opt.frame_range = [90:150];
% get two choice decoders for go trials and no-go trials separately
choice_decods = struct();

dpca_fd_names ={ {'stim_1_correct','stim_1_incorrect','stim_2_incorrect','stim_2_correct'},...
    {'stim_1_correct','stim_1_incorrect','stim_2_incorrect','stim_2_correct'}};  

decod_fd_names ={ {'stim_1_correct','stim_1_incorrect'},...
    {'stim_2_incorrect','stim_2_correct'}};  % first field will be positive

decod_cell_types = {'choice'}; 
% 
for d = 1:numel(decod_fd_names)
    % dpca
     if d<=numel(decod_cell_types) % not going to compute again if using same cells
        dpca_opt.cell_idx = cell_idx_struct.(decod_cell_types{d});    
        dpca_opt.fd_names = dpca_fd_names{d};
        
        
        switch opt.method
            case 'dpca'
                [dpca_struct,traces_idx_struct] = get_dpca_traj_brief(baseline_train_cell_struct,dpca_opt.cell_idx,dpca_opt.fd_names,dpca_opt,'frame_range',dpca_opt.frame_range);
                dpca_struct.std = ones(size(dpca_struct.mean));
                traj_struct = dpca_struct.traj_struct;
%                 pcs_to_use = unique(cell2mat(struct2cell(dpca_struct.PCidx)')); % get condition dependent PCs
                pcs_to_use = find(dpca_struct.whichMarg~=3);
pcs_to_use =  find(dpca_struct.whichMarg==1);
            case 'pca'
                % % pca
                [traj_struct,dpca_struct,traj_trial_idx] = get_cell_pca(baseline_train_cell_struct(dpca_opt.cell_idx),dpca_opt);
                dpca_struct.W = dpca_struct.U(:,1:5); % usinng 'U' from svd as weights
                pcs_to_use = 1:3;
                
            case 'cd'
                % coding direction
                choice_fds = {'stim_1_correct','stim_2_correct'};
                cd_opt = opt;
                cd_opt.trial_color = trial_color;
                coding_direct = cell2mat(arrayfun(@(x)nanmean(mean(baseline_train_cell_struct(x).(choice_fds{1})(opt.sta_avg_frames,:),2)),dpca_opt.cell_idx,'UniformOutput', false))...
                    - cell2mat(arrayfun(@(x)nanmean(mean(baseline_train_cell_struct(x).(choice_fds{2})(opt.sta_avg_frames,:),2)),dpca_opt.cell_idx,'UniformOutput', false));
                % project to coding direction
                traj_struct = get_cd_projection(baseline_train_cell_struct(dpca_opt.cell_idx),coding_direct,dpca_opt.fd_names );
                pcs_to_use = 1;
                dpca_struct.W = ones(numel(dpca_opt.cell_idx),1);
                
                
            case 'diff'
                catfun = @(C)mean(reshape(cell2mat(C),size(C{1},1),size(C{1},2),size(C,2)),3);
                % difference between texture-coding ensembles
                for fd = dpca_opt.fd_names
                    this_num_trials = size(baseline_train_cell_struct(1).(fd{:}),2);
                    
                    for t = 1:this_num_trials
                        traj_struct.(fd{:})= catfun({baseline_train_cell_struct(cell_idx_struct.('whisk_tex1')).(fd{:})})'-...
                            catfun({baseline_train_cell_struct(cell_idx_struct.('whisk_tex2')).(fd{:})})';
                    end
                end
                
                dpca_struct.W = [ones(numel(cell_idx_struct.('whisk_tex1')),1); -ones(numel(cell_idx_struct.('whisk_tex2')),1)];
                dpca_struct.cell_idx = dpca_opt.cell_idx;
                dpca_struct.mean = zeros(1,numel(dpca_struct.cell_idx));
                dpca_struct.std = ones(1,numel(dpca_struct.cell_idx));
                pcs_to_use = 1;
        end
        num_pcs = numel(pcs_to_use);
        
        plot_pop_vectors(traj_struct,dpca_opt.fd_names,num_pcs,dpca_opt,...
            'plot_ylabel','PC level','plot_num_cols',2);
        
        
        for fd = dpca_opt.fd_names
            this_traj_traces = traj_struct.(fd{:})(:,:,pcs_to_use);
            %     this_traces = nan(size(this_traj_traces,1),size(this_traj_traces,2),num_pcs+1);% append neuropil trace to trajectory background_sta_traces
            %     this_traces(:,:,1:num_pcs) = this_traj_traces;
            %     this_traces(:,:,num_pcs+1)= baseline_train_bg_struct.(fd{:});
            traj_struct.(strrep(fd{:},'_nonphoto',''))=this_traj_traces;
        end
        
        %     traj_struct.('go_trials') = cat(1,traj_struct.stim_1_incorrect,traj_struct.stim_2_correct);
        %     traj_struct.('nogo_trials') = cat(1,traj_struct.stim_1_correct,traj_struct.stim_2_incorrect);
     end
    
    
    
choice_opt.fd_names = decod_fd_names{d}; %correct choice will be positive
choice_opt.frame_range = 1:opt.trial_length;

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
% export_fig([fig_save_path filesep 'ChoiceDecoderPerform_' strrep(caiman_file,'.mat','.png')])
choice_weights = dpca_struct.W(:,pcs_to_use)* choice_struct.B(:,2:end)';
choice_thresh = -choice_struct.B(:,1);
[choice_norm_weights,choice_norm_thresh] = get_norm_weights(choice_weights,choice_thresh,dpca_struct.mean,dpca_struct.std); % called in CompareDecoderWithAnimal

choice_decods(d).choice_struct = choice_struct;
choice_decods(d).choice_norm_weights = choice_norm_weights;
choice_decods(d).choice_norm_thresh = choice_norm_thresh;
choice_decods(d).fd_names = choice_opt.fd_names ;
choice_decods(d).cell_idx = dpca_opt.cell_idx;
end

% project to decoders
decod_proj_struct = struct();
plot_fds = {'stim_1_correct','stim_1_incorrect','stim_2_correct','stim_2_incorrect'};
for f = 1:numel(plot_fds)
    fd = plot_fds{f};
    this_decod = strsplit(fd,'_');
    this_decod = str2double(this_decod{2});
    norm_thresh = choice_decods(this_decod).choice_norm_thresh;
    norm_weights = zeros(num_cells,1);
    norm_weights(choice_decods(this_decod).cell_idx ) = choice_decods(this_decod).choice_norm_weights;
    [decod_proj_struct] = get_projections(baseline_train_cell_struct,norm_weights,{fd},'bias',-norm_thresh,'IS_CELL_STRUCT',1,'proj_struct',decod_proj_struct);
end
%% plot decoders performance
for d = 1:numel(decod_fd_names)
    figure
    plot_binary_decoder(choice_decods(d).choice_struct,choice_opt)
    suptitle(['Tex ' num2str(d) ' decoder'])
end
%% 4. get texture decoder
tex_opt = dpca_opt;
tex_opt.frames_to_avg = [130:145];% before motor stops
tex_opt.frames_to_train = round([1:1:150]/tex_opt.bin_size);
tex_opt.min_frames = 10;
tex_opt.Nstd = 1.5;
tex_opt.IF_FRAMEWISE = 0;
tex_opt.Fs = 30;
tex_opt.fd_names = {'stim_1_correct','stim_2_correct'};
tex_opt.frame_range = 1:opt.trial_length;
tex_struct = {};
tex_proj_struct = {};
disp('Rnning texture decoder...')
tic
tex_struct =  get_binary_classifier(tex_struct,traj_struct, tex_opt,...
    'IF_CROSSVAL',1,'IF_FRAMEWISE',tex_opt.IF_FRAMEWISE,'fd_names',tex_opt.fd_names);
[tex_proj_struct] = get_projections(traj_struct,tex_struct.B(:,2:end)',dpca_opt.fd_names,'proj_struct',tex_proj_struct,'bias',tex_struct.B(:,1));
plot_pop_vectors(tex_proj_struct,dpca_opt.fd_names,1,tex_opt,...
        'plot_ylabel','Texture decoder projection')
[ tex_struct ] =  get_binary_decoder_disc_time( tex_proj_struct, tex_struct,...
    tex_opt.fd_names,tex_opt,'IF_FRAMEWISE',tex_opt.IF_FRAMEWISE,'threshold',0);
toc
disp('Done')

figure; 
hold on
plot_binary_decoder(tex_struct,tex_opt)
suptitle('Texture decoder')
% export_fig([fig_save_path filesep 'ChoiceDecoderPerform_' strrep(caiman_file,'.mat','.png')])
tex_weights = dpca_struct.W(:,pcs_to_use)* tex_struct.B(:,2:end)';
tex_thresh = -tex_struct.B(:,1);
[tex_norm_weights,tex_norm_thresh] = get_norm_weights(tex_weights,tex_thresh,dpca_struct.mean,dpca_struct.std); % called in CompareDecoderWithAnimal

%% distance from correct trajectory (noisy)
[tex_proj_struct] = get_projections(baseline_train_cell_struct(tex_opt.cell_idx),tex_norm_weights,plot_fds,'bias',-tex_norm_thresh,'IS_CELL_STRUCT',1);
plot_pop_vectors(tex_proj_struct,plot_fds,1,tex_opt,...
        'plot_ylabel','Tex. decoder projection')

ref_fds = {'stim_1_correct','stim_2_correct'};
[rel_dist_traces] = get_dist_to_ref(tex_proj_struct,ref_fds,plot_fds);
    
figure('name','selectivity: relative distance (baseline session)','position',[100 100 800 400])
values = structfun(@(x)mean(x(:,tex_opt.frames_to_avg),2),rel_dist_traces,'UniformOutput',false);
fd_colors =  cell2mat(cellfun(@(f)getfield(trial_color,f),fields(values),'UniformOutput',false));
scatter_cmp_conditions(values,[],...
    1,fd_colors,'plot_stats',1,'connect_scatter',0,'BriefXlabel',0,'ShowMeanInXlabel',1,'add_jitter',1,'xtick_label',plot_fds_names);
axis square
xtickangle(45)
ylabel('Rel. distance')

%% =============Simulation of error detection using online decoder  =================
%% load config file: baseline output struct (file named as 'OutputParams'
% reproduce online decoder accuracy estimated from baseline session
% skip this for post hoc simulation

[baseline_file,baseline_path] = uigetfile('*.mat','Select baseline OutputParams',caiman_path);
baseline_output = load(fullfile(baseline_path,baseline_file));
baseline_output = baseline_output.output;
disp(['Loaded file :',fullfile(baseline_path,baseline_file)])
% trajectory computed using raw_sta_traces (reproduce online traj)
% weights from online output
norm_thresh = baseline_output.trigger_thresh;
norm_weights = cell2mat({cell_struct(:).weight})';

plot_fds =  {'stim_1_correct','stim_1_incorrect','stim_2_correct','stim_2_incorrect'};
plot_fds_names = {'CR','FA','Hit','Miss'};
choice_opt.frames_to_avg = 90:120;
[decod_proj_struct] = get_projections(baseline_raw_cell_struct,norm_weights,plot_fds,'bias',-norm_thresh,'IS_CELL_STRUCT',1,'xtick_label',plot_fds_names);
%% 1. baseline projection traces 
plot_fds = {'stim_1_correct','stim_1_incorrect','stim_2_correct','stim_2_incorrect'};% {'stim_1_correct','stim_1_incorrect','stim_2_correct','stim_2_incorrect'};
plot_fds_names = {'CR','FA','Hit','Miss'};
plot_pop_vectors(decod_proj_struct,plot_fds,1,opt,...
    'plot_ylabel','Projection','plot_num_cols',4,'IF_PLOT_RAW_ONLY',0)
% xlim([3 6])
ylabel('Decoder projection')
figure('name','choice decoder projection (baseline session)','position',[100 100 800 400])
subplot(1,3,1)
values = structfun(@(x)mean(x(:,choice_opt.frames_to_avg),2),decod_proj_struct,'UniformOutput',false);
fd_colors =  cell2mat(cellfun(@(f)getfield(trial_color,f),fields(values),'UniformOutput',false));
scatter_cmp_conditions(values,[],...
    1,fd_colors,'plot_stats',1,'connect_scatter',0,'BriefXlabel',0,'ShowMeanInXlabel',1,'add_jitter',1,'xtick_label',plot_fds_names);
axis square
xtickangle(45)
ylabel('Decoder projection')
title('baseline session')

% expected hit and fa rate of error detection
decod_accuracy = struct();
decod_accuracy.CR = numel(find(values.stim_1_correct>0))./numel(values.stim_1_correct);
decod_accuracy.FA = numel(find(values.stim_1_incorrect<0))./numel(values.stim_1_incorrect);
decod_accuracy.Hit = numel(find(values.stim_2_correct<0))./numel(values.stim_2_correct);
decod_accuracy.Miss = numel(find(values.stim_2_incorrect>0))./numel(values.stim_2_incorrect);

overall = mean(cell2mat(struct2cell(decod_accuracy)));
subplot(1,3,2)
scatter_cmp_conditions(decod_accuracy,[],...
    1,fd_colors,'connect_scatter',0,'BriefXlabel',1,'ShowMeanInXlabel',0,'add_jitter',0,'xtick_label',plot_fds_names);
axis square
xtickangle(45)
ylabel('Decoder accuracy')
title(['overall accuracy: ' num2str(overall*100,'%.2f') '%'])

% export_fig([fig_save_path filesep 'BaselineTwoChoiceDecodersProjWhiskingTrials' '.png'])

% run decoder on baseline session traces to get error trials (simulation)
plot_trial_indices = baseline_trial_indices;
% monitor_frames = opt.sta_baseline_frames+baseline_output.trigger_frames;
monitor_frames =  opt.sta_baseline_frames + [100 120];
[decod_accuracy,overall] = get_expected_online_decod_accuracy(decod_proj_struct,plot_trial_indices,monitor_frames);
% figure('name','expected decoder accuracy','position',[100 100 800 400])
subplot(1,3,3)
title(['overall accuracy: ' num2str(overall*100,'%.2f') '%'])
fd_colors =  cell2mat(cellfun(@(f)getfield(trial_color,f),plot_fds','UniformOutput',false));
scatter_cmp_conditions(decod_accuracy,[],...
    1,fd_colors,'connect_scatter',0,'BriefXlabel',1,'ShowMeanInXlabel',0,'add_jitter',0);
ylim([0 1])
axis square
xtickangle(45)
ylabel('Expected online decoder accuracy')

suptitle('Choice decoder projections (baseline session)')

% save([session_path '\analysis_files\baseline_dummy_indices'],'dummy_indices');
%% plot dummy stim trajs
trial_types = fields(dummy_indices);
for i = 1:numel(trial_types)
    this_fd = trial_types{i};
    this_idx = dummy_indices.(this_fd);
    this_idx = this_idx(this_idx<=min([tot_num_trials,num_sta_traces]));
%     this_idx = setdiff(this_idx,easy_trial_idx);
    if ~isempty(this_idx)
        for c = 1:num_cells
            baseline_whisk_cell_struct(c).(this_fd) = baseline_cell_struct(c).raw_sta_traces( this_idx,:)';
        end
    end
end

% [decod_proj_struct] = get_projections(baseline_raw_cell_struct,norm_weights,fields(dummy_indices),'bias',-norm_thresh,'IS_CELL_STRUCT',1);
% project to decoders
trial_types = {{'stim_1_dummyphoto_correct','stim_1_dummyphoto_incorrect'},...
    {'stim_2_dummyphoto_correct','stim_2_dummyphoto_incorrect'}};

decod_proj_struct = struct();
for d = 1:numel(decod_fd_names)
norm_thresh = choice_decods(d).choice_norm_thresh;
norm_weights = zeros(num_cells,1);
norm_weights(cell_idx_struct.(opt.trigger_idx_fd)) = choice_decods(d).choice_norm_weights;
[decod_proj_struct] = get_projections(baseline_whisk_cell_struct,norm_weights,trial_types{d},...
    'bias',-norm_thresh,'IS_CELL_STRUCT',1,'proj_struct',decod_proj_struct);
end


plot_pop_vectors(decod_proj_struct,fields(dummy_indices),1,opt,...
    'plot_ylabel','Projection','plot_num_cols',4,'IF_PLOT_RAW_ONLY',1)

%% 2. test online decoder performance using condition session (simulation)
% pool all trials, test on frames before photostim enabled
trial_types =  {'stim_1_correct','stim_1_incorrect','stim_2_correct','stim_2_incorrect'};

whisk_trial_indices.stim_1_correct = unique([whisk_trial_indices.stim_1_photostim_correct whisk_trial_indices.stim_1_nonphotostim_correct whisk_trial_indices.stim_1_dummyphotostim_correct]);
whisk_trial_indices.stim_1_incorrect = unique([whisk_trial_indices.stim_1_photostim_incorrect whisk_trial_indices.stim_1_nonphotostim_incorrect whisk_trial_indices.stim_1_dummyphotostim_incorrect]);
whisk_trial_indices.stim_2_correct = unique([whisk_trial_indices.stim_2_photostim_correct whisk_trial_indices.stim_2_nonphotostim_correct whisk_trial_indices.stim_2_dummyphotostim_correct]);
whisk_trial_indices.stim_2_incorrect = unique([whisk_trial_indices.stim_2_photostim_incorrect whisk_trial_indices.stim_2_nonphotostim_incorrect whisk_trial_indices.stim_2_dummyphotostim_incorrect]);

% % not using trials with photostim
whisk_trial_indices.stim_1_correct = unique([ whisk_trial_indices.stim_1_nonphotostim_correct whisk_trial_indices.stim_1_dummyphotostim_correct]);
whisk_trial_indices.stim_1_incorrect = unique([whisk_trial_indices.stim_1_nonphotostim_incorrect whisk_trial_indices.stim_1_dummyphotostim_incorrect]);
whisk_trial_indices.stim_2_correct = unique([whisk_trial_indices.stim_2_nonphotostim_correct whisk_trial_indices.stim_2_dummyphotostim_correct]);
whisk_trial_indices.stim_2_incorrect = unique([ whisk_trial_indices.stim_2_nonphotostim_incorrect whisk_trial_indices.stim_2_dummyphotostim_incorrect]);



trial_indices.stim_1_correct = unique([trial_indices.stim_1_nonphotostim_correct trial_indices.stim_1_dummyphotostim_correct]);
trial_indices.stim_1_incorrect = unique([trial_indices.stim_1_photostim_incorrect trial_indices.stim_1_nonphotostim_incorrect trial_indices.stim_1_dummyphotostim_incorrect]);
trial_indices.stim_2_correct = unique([trial_indices.stim_2_photostim_correct trial_indices.stim_2_nonphotostim_correct trial_indices.stim_2_dummyphotostim_correct]);
trial_indices.stim_2_incorrect = unique([trial_indices.stim_2_photostim_incorrect trial_indices.stim_2_nonphotostim_incorrect trial_indices.stim_2_dummyphotostim_incorrect]);

condition_whisk_cell_struct = struct();
for i = 1:numel(trial_types)
    this_fd = trial_types{i};
    this_idx = whisk_trial_indices.(this_fd);
    if ~isempty(this_idx)
        for c = 1:num_cells
            condition_whisk_cell_struct(c).(this_fd) = condition_cell_struct(c).raw_sta_traces( this_idx,:)';
        end
    end
    condition_whisk_trial_indices.(strrep(this_fd,'_nonphoto','')) = this_idx;
end
plot_fds = trial_types;
condition_decod_proj_struct = struct();
for f = 1:numel(trial_types)
    fd = plot_fds{f};
    this_decod = strsplit(fd,'_');
    this_decod = str2double(this_decod{2});
    norm_thresh = choice_decods(this_decod).choice_norm_thresh;
    norm_weights = zeros(num_cells,1);
    norm_weights(choice_decods(this_decod).cell_idx) = choice_decods(this_decod).choice_norm_weights;
    [condition_decod_proj_struct] = get_projections(condition_whisk_cell_struct,norm_weights,{fd},'bias',-norm_thresh,'IS_CELL_STRUCT',1,'proj_struct',condition_decod_proj_struct);
end

%% check dpca traces in condition session
[condition_traj_struct] = get_traj_struct(condition_whisk_cell_struct(dpca_opt.cell_idx),dpca_struct.W(:,pcs_to_use),plot_fds,'cell_mean',[]);
plot_pop_vectors(condition_traj_struct,plot_fds,num_pcs,tex_opt,...
    'plot_ylabel','PC' ,'IF_PLOT_AVG_ONLY',1)

plot_pop_vectors(traj_struct,plot_fds,num_pcs,tex_opt,...
    'plot_ylabel','PC' ,'IF_PLOT_AVG_ONLY',1)

condition_choice_proj_struct = struct();
for f = 1:numel(plot_fds)
    fd = plot_fds{f};
    this_decod = strsplit(fd,'_');
    this_decod = str2double(this_decod{2});     
    [condition_choice_proj_struct] = get_projections(condition_traj_struct,choice_decods(this_decod).choice_struct.B(:,2:end)',{fd},'proj_struct',condition_choice_proj_struct,'bias',choice_decods(this_decod).choice_struct.B(:,1));
end

condition_decod_proj_struct = condition_choice_proj_struct;
%%
% figure;
% condition_decod_proj_struct = condition_choice_proj_struct;
plot_pop_vectors(condition_decod_proj_struct,plot_fds,1,opt,...
        'plot_ylabel','Projection','plot_num_cols',4,'IF_PLOT_RAW_ONLY',0)

%% plot decoder accuracy (condition)
figure('name','choice decoder projection (condition session)','position',[100 100 800 400])
subplot(1,3,1)

values = structfun(@(x)mean(x(:,choice_opt.frames_to_avg),2),condition_decod_proj_struct,'UniformOutput',false);
fd_colors =  cell2mat(cellfun(@(f)getfield(trial_color,strrep(f,'_var_5','')),fields(values),'UniformOutput',false));
scatter_cmp_conditions(values,[],...
    1,fd_colors,'plot_stats',1,'connect_scatter',0,'BriefXlabel',0,'ShowMeanInXlabel',1,'add_jitter',1,'xtick_label',{'CR','FA','Hit','Miss'});
axis square
xtickangle(45)
ylabel('Decoder projections')
title('condition session')


% expected hit and fa rate of error detection
subplot(1,3,2)
decod_accuracy = struct();
decod_accuracy.CR = numel(find(values.stim_1_correct>0))./numel(values.stim_1_correct);
decod_accuracy.FA = numel(find(values.stim_1_incorrect<0))./numel(values.stim_1_incorrect);
decod_accuracy.Hit = numel(find(values.stim_2_correct<0))./numel(values.stim_2_correct);
decod_accuracy.Miss = numel(find(values.stim_2_incorrect>0))./numel(values.stim_2_incorrect);

overall = mean(cell2mat(struct2cell(decod_accuracy)));
scatter_cmp_conditions(decod_accuracy,[],...
    1,fd_colors,'connect_scatter',0,'BriefXlabel',0,'ShowMeanInXlabel',1,'add_jitter',0,'xtick_label',plot_fds_names);
axis square
xtickangle(45)
ylabel('Decoder accuracy')
title(['overall accuracy: ' num2str(overall*100,'%.2f') '%'])



[decod_accuracy,overall] = get_expected_online_decod_accuracy(condition_decod_proj_struct,whisk_trial_indices,monitor_frames);
subplot(1,3,3)
title(['overall accuracy: ' num2str(overall*100,'%.2f') '%'])
fd_colors =  cell2mat(cellfun(@(f)getfield(trial_color,f),plot_fds','UniformOutput',false));
scatter_cmp_conditions(decod_accuracy,[],...
    1,fd_colors,'connect_scatter',0,'BriefXlabel',1,'ShowMeanInXlabel',0,'add_jitter',0);
ylim([0 1])
axis square
xtickangle(45)
ylabel('Expected online decoder accuracy')
suptitle('Choice decoder projections (condition session)')

%% project to texture decoder 
[condition_tex_proj_struct] = get_projections(condition_whisk_cell_struct(tex_opt.cell_idx),tex_norm_weights,plot_fds,'bias',-tex_norm_thresh,'IS_CELL_STRUCT',1);
plot_pop_vectors(condition_tex_proj_struct,plot_fds,1,tex_opt,...
    'plot_ylabel','Tex. decoder projection','IF_PLOT_AVG_ONLY',1)
suptitle('condition session')
plot_pop_vectors(tex_proj_struct,plot_fds,1,tex_opt,...
    'plot_ylabel','Tex. decoder projection','IF_PLOT_AVG_ONLY',1)
suptitle('Baseline session')


%% test texture decoder using condition session
[condition_decod_proj_struct] = get_projections(condition_whisk_cell_struct(tex_opt.cell_idx),tex_norm_weights,plot_fds,'bias',-tex_norm_thresh,'IS_CELL_STRUCT',1);
plot_pop_vectors(condition_decod_proj_struct,plot_fds,1,tex_opt,...
        'plot_ylabel','Tex. decoder projection')

ref_fds = {'stim_1_correct','stim_2_correct'};
[rel_dist_traces] = get_dist_to_ref(condition_decod_proj_struct,ref_fds,plot_fds);
    
figure('name','selectivity: relative distance (condition session)','position',[100 100 800 400])
values = structfun(@(x)mean(x(:,tex_opt.frames_to_avg),2),rel_dist_traces,'UniformOutput',false);
fd_colors =  cell2mat(cellfun(@(f)getfield(trial_color,f),fields(values),'UniformOutput',false));
scatter_cmp_conditions(values,[],...
    1,fd_colors,'plot_stats',1,'connect_scatter',0,'BriefXlabel',0,'ShowMeanInXlabel',1,'add_jitter',1);
axis square
xtickangle(45)
ylabel('Rel. distance')


plot_pop_vectors(condition_decod_proj_struct,plot_fds,1,opt,...
        'plot_ylabel','Projection','plot_num_cols',4,'IF_PLOT_RAW_ONLY',0)


%% 3. test online decoder performance using control session (simulation)
% get non-photo trials from control session (var 5)
trial_types = fields(ctr_whisk_trial_indices);
ctr_whisk_cell_struct = struct();
for i = 1:numel(trial_types)
    this_fd = trial_types{i};
    this_idx = ctr_whisk_trial_indices.(this_fd);
    if ~isempty(this_idx)
        for c = 1:num_cells
            ctr_whisk_cell_struct(c).(strrep(this_fd,'_nonphoto','')) = ctr_cell_struct(c).raw_sta_traces( this_idx,:)';
        end
    end
    ctr_whisk_trial_indices.(strrep(this_fd,'_nonphoto','')) = this_idx;
end

ctr_plot_fds = {'stim_1_var_5_correct','stim_1_var_5_incorrect','stim_2_var_5_correct','stim_2_var_5_incorrect'};
ctr_decod_proj_struct = struct();
for f = 1:numel(ctr_plot_fds)
    fd = ctr_plot_fds{f};
    this_decod = strsplit(fd,'_');
    this_decod = str2double(this_decod{2});
    norm_thresh = choice_decods(this_decod).choice_norm_thresh;
    norm_weights = zeros(num_cells,1);
    norm_weights(dpca_opt.cell_idx) = choice_decods(this_decod).choice_norm_weights;
    [ctr_decod_proj_struct] = get_projections(ctr_whisk_cell_struct,norm_weights,{fd},'bias',-norm_thresh,'IS_CELL_STRUCT',1,'proj_struct',ctr_decod_proj_struct);
end

plot_pop_vectors(ctr_decod_proj_struct,ctr_plot_fds,1,opt,...
        'plot_ylabel','Projection','plot_num_cols',4,'IF_PLOT_RAW_ONLY',1)

% figure;
% plot(ctr_raw_cell_struct(1).stim_1_var_5_correct)
% suptitle('Decoder projections')
figure('name','choice decoder projection (baseline session)','position',[100 100 800 400])
values = structfun(@(x)mean(x(:,opt.sta_avg_frames),2),ctr_decod_proj_struct,'UniformOutput',false);
fd_colors =  cell2mat(cellfun(@(f)getfield(trial_color,strrep(f,'_var_5','')),fields(values),'UniformOutput',false));
scatter_cmp_conditions(values,[],...
    1,fd_colors,'connect_scatter',0,'BriefXlabel',0,'ShowMeanInXlabel',1,'add_jitter',1,'xtick_label',{'CR','FA','Hit','Miss'});
axis square

suptitle('Choice decoder projections (control session)')


%% ===================    CONDITION SESSION   reproduce  ===========================
%% plot photostim ensemble texture sta traces
cell_types = fields(cell_idx_struct);
plot_fds = {'stim_1_correct','stim_1_incorrect','stim_2_correct','stim_2_incorrect'};
plot_fds_names = {'CR','FA','Hit','Miss'};
for f = 1:numel(cell_types)
    for ff = 1:numel(plot_fds)
        fd = plot_fds{ff};
        ensemble_mean.(cell_types{f}).(fd) = [];
        ensemble_traces.(cell_types{f}).(fd) = [];
    end
end

[ensemble_traces,ensemble_mean,ensemble_traces_raw,ensemble_mean_raw] = get_ensemble_mean(ensemble_mean,ensemble_traces,condition_whisk_cell_struct,cell_idx_struct,plot_fds,opt);

figure('name','ensemble traces','units','normalized','outerposition',[0 0 1 .8])
plot_cell_types = {'all','stim1','stim2','choice1','choice2'};

trace_plot_fds = {'stim_1_correct','stim_2_correct','stim_1_incorrect','stim_2_incorrect'};
trace_names = {'CR','Hit','FA','Miss'};

this_xticks = ([1:opt.trial_length]-opt.sta_gocue_frame)./opt.frame_rate;
plot_count = 1;
num_plots = numel(plot_cell_types);
for p = 1:numel(plot_cell_types)
    subplot(2,num_plots,plot_count)
    
    hold on
    
    cell_id = plot_cell_types{p};
    this_num_cells = numel(cell_idx_struct.(cell_id));
    for f = 1:numel(trace_plot_fds)
        fd = trace_plot_fds{f};
        this_traces = ensemble_traces_raw.(cell_id).(fd)';
        shadedErrorBar(this_xticks,mean(this_traces,1), std(this_traces,[],1),{'color',trial_color.(fd),'linewidth',2},0.1);
    end
    ylabel([ strrep(cell_id,'_','-') ' ens. avg., ' num2str(this_num_cells) ' cells' ])
    xlabel('Time from go-cue (sec)')
    ylim([-0.5,5])
    xlim([-4,4])

    plot([0 0],ylim,':','color',[.5 .5 .5],'linewidth',2) % go-cue
    plot([-1 0],[0 0],'color','black','linewidth',2) % sta average window (withold, 1sec before gocue)
    axis square
    
    subplot(2,num_plots,plot_count+num_plots)

    this_values = ensemble_mean_raw.(cell_id);
    fd_colors =  cell2mat(cellfun(@(f)getfield(trial_color,f),fields(this_values),'UniformOutput',false));
    
    scatter_cmp_conditions(this_values,[],...
        1,fd_colors,'connect_scatter',0,'BriefXlabel',1,'ShowMeanInXlabel',1,'add_jitter',1,...
        'plot_stats',1,'test_type','ranksum','xtick_label',plot_fds_names);
    ylabel([ strrep(cell_id,'_','-') ' cell avg.' ])
    axis square
     plot_count =plot_count+ 1;
   
end
% export_fig([fig_save_path filesep 'OnlineTexCellsEnsembleMean' '.png'])

%% Plot cells of interest on FOV
figure('name','trigger targets check','position',[100 100 2400 800]);
for dummy = 1
ax = subtightplot(1,4,1);
plot_value_in_rois( cell_struct, 'cnn_prediction',[256 256],ax,...
    'IF_NORM_PIX',1,'IF_CONTOUR',0,'colorlut',flipud(gray));
set(gca,'Ydir','reverse')
title('CNN prediction score')

% mark texture cells on FOV
% mark_cells_on_fov(cell_struct,cell_idx_struct.(opt.trigger_idx_fd),[.5 .5 .5],'MarkerSize',300,'Linewidth',2)
% unique_target_fds = unique(opt.target_idx_fd);
% for f = 1:numel(unique_target_fds)
%     fd = strsplit(unique_target_fds{f},'_');
%     fdd = fd{end};
%     mark_cells_on_fov(cell_struct,cell_idx_struct.(unique_target_fds{f}),trial_color.(fdd),'MarkerSize',500,'Linewidth',2)
% end

ax = subtightplot(1,4,2);
plot_value_in_rois( cell_struct, 'correct_stimAUC_zscore',[256 256],ax,...
    'IF_NORM_PIX',1,'IF_CONTOUR',0,'zlimit',[-3 3],'colorlut',flipud(b2r(-3 ,3)));
set(gca,'Ydir','reverse')
title('Texture selectivity (auc zscore)')

% mark texture cells on FOV
% mark_cells_on_fov(cell_struct,cell_idx_struct.(opt.trigger_idx_fd),[.5 .5 .5],'MarkerSize',300,'Linewidth',2)
unique_target_fds = {'tex1','tex2'};
for f = 1:numel(unique_target_fds)
    fd = strsplit(unique_target_fds{f},'_');
    fdd = fd{end};
    mark_cells_on_fov(cell_struct,cell_idx_struct.(unique_target_fds{f}),trial_color.(fdd),'MarkerSize',500,'Linewidth',2)
end

ax = subtightplot(1,4,3);
plot_value_in_rois( cell_struct, 'weight',[256 256],ax,...
    'IF_NORM_PIX',1,'IF_CONTOUR',0);
set(gca,'Ydir','reverse')
title('Decoder weight')

% mark texture cells on FOV
mark_cells_on_fov(cell_struct,cell_idx_struct.(opt.trigger_idx_fd),[.5 .5 .5],'MarkerSize',300,'Linewidth',2);


ax = subtightplot(1,4,4);
plot_value_in_rois( cell_struct, 'photo_auc_zscore',[256 256],ax,...
    'IF_NORM_PIX',1,'IF_CONTOUR',0,'zlimit',[-4 4]);
set(gca,'Ydir','reverse')
title('Photo response (auc zscore)')

% mark texture cells on FOV
% mark_cells_on_fov(cell_struct,cell_idx_struct.(opt.trigger_idx_fd),[.5 .5 .5],'MarkerSize',300,'Linewidth',2)
unique_target_fds = unique(opt.target_idx_fd);
for f = 1:numel(unique_target_fds)
    fd = strsplit(unique_target_fds{f},'_');
    fdd = fd{end};
    mark_cells_on_fov(cell_struct,cell_idx_struct.(unique_target_fds{f}),trial_color.(fdd),'MarkerSize',500,'Linewidth',2)
end
end
export_fig([fig_save_path filesep 'OnlineIdentifiedCells' '.pdf'])

%% mask texture selectivity on fov
figure
ax = subplot(1,2,1);
plot_cell_idx = tex_output_struct.cell_idx_struct.cnn_above_thresh;
plot_cell_idx = 1:size(cell_struct,2);

fov_img = plot_value_in_rois( cell_struct(plot_cell_idx), [],[256 256],ax,...
    'IF_NORM_PIX',1,'IF_CONTOUR',0,'colorlut',gray);
set(gca,'Ydir','reverse')
colormap(gray)

tex_img = plot_value_in_rois( cell_struct(plot_cell_idx), 'correct_stimAUC_zscore',[256 256],ax,...
    'IF_NORM_PIX',0,'IF_CONTOUR',0,'zlimit',[-4 4],'colorlut',flipud(b2r(-4 ,4)));
set(gca,'Ydir','reverse')



texlabel_img = zeros(size(tex_img));
texlabel_img(tex_img<-opt.N) = -1;
texlabel_img(tex_img>opt.N) = 1;

close

figure('position',[100 100 1000 1000])
subplot(1,2,1)
hold on
C = labeloverlay(fov_img,categorical(texlabel_img),...
    'Colormap',[1 0 0;0 0 0;0 0 1],'Transparency',0.7);
imshow(C)
% mark targeted cells on FOV
unique_target_fds = unique(opt.target_idx_fd);
for f = 1:numel(unique_target_fds)
    fd = strsplit(unique_target_fds{f},'_');
    fdd = fd{end};
    hold on
    mark_cells_on_fov(cell_struct,intersect(plot_cell_idx,cell_idx_struct.(unique_target_fds{f})),trial_color.(fdd),'MarkerSize',200,'Linewidth',2)
end


ax = subplot(1,2,2);
hold on
plot_value_in_rois( cell_struct(plot_cell_idx), 'photo_auc_zscore',[256 256],ax,...
    'IF_NORM_PIX',1,'IF_CONTOUR',0,'zlimit',[-2 2],'colorlut',gray);
set(gca,'Ydir','reverse')

% mark targeted cells on FOV
unique_target_fds = unique(opt.target_idx_fd);
for f = 1:numel(unique_target_fds)
    fd = strsplit(unique_target_fds{f},'_');
    fdd = fd{end};
    hold on
    mark_cells_on_fov(cell_struct,intersect(plot_cell_idx,cell_idx_struct.(unique_target_fds{f})),trial_color.(fdd),'MarkerSize',200,'Linewidth',2)
end

unique_target_fds = {'tex1','tex2'};
for f = 1:numel(unique_target_fds)
    fd = strsplit(unique_target_fds{f},'_');
    fdd = fd{end};
    hold on
    mark_cells_on_fov(cell_struct,intersect(plot_cell_idx,cell_idx_struct.(unique_target_fds{f})),trial_color.(fdd),'MarkerSize',150,'Linewidth',1)
end

% export_fig([fig_save_path filesep 'OnlineIdentifiedCells' '.pdf'])
%% get online trajectory and weights
online_weights = condition_caiman_data.ROIw;
online_thresh = condition_caiman_data.ROIsumThresh;
online_sd = condition_caiman_data.sd_level;
online_bs = condition_caiman_data.bs_level;
online_w = condition_caiman_data.ROIw;
disp('got online trajectory and weights')
%% STA of online trajectory
online_traj = condition_caiman_data.online_traj;
online_traj = online_traj./max(abs(online_traj));
num_frames = condition_caiman_data.t_cnm;
skip_frames = condition_caiman_data.frames_skipped + condition_caiman_data.t_init+1;
tot_frames = num_frames + numel(skip_frames);
caiman_frames = setdiff([1:tot_frames],skip_frames);

traj_struct = struct();
temp_trace = nan(1,tot_frames);
temp_trace(caiman_frames) = online_traj(1,1:num_frames);
%  online_traj = fillmissing(temp_trace,'linear');
online_traj = temp_trace;
sens_stim_frames = condition_caiman_data.sensory_stim_frames+condition_caiman_data.t_init;
this_sens_stim_frames = sens_stim_frames+opt.gocue_frame;
[~,~,~,traj_struct.sta_traces,~,~,traj_struct.sta_trace] =...
    make_sta_from_traces(online_traj,this_sens_stim_frames,opt.sta_pre_frames,opt.sta_post_frames,1:opt.sta_baseline_frames);

% Sort STAs by different trial types
trial_types = fields(trial_indices);
num_sta_traces = size(traj_struct(1).sta_traces,1);

for i = 1:numel(trial_types)
    this_fd = trial_types{i};
%     if contains(this_fd,'stim_1_')
%         this_offset = -online_thresh;
%     elseif contains(this_fd,'stim_2_')
%         this_offset = online_thresh;   
%     end
    this_idx = trial_indices.(this_fd);
    this_idx = this_idx(this_idx<=min([tot_num_trials,num_sta_traces]));
%     this_idx = setdiff(this_idx,easy_trial_idx);
    if ~isempty(this_idx)
        for c = 1:num_cells
            traj_struct.(this_fd) = traj_struct.sta_traces(this_idx,:);
        end
    end
end

disp('sorted traj sta')

%% add condition session traces to cell struct
num_frames = condition_caiman_data.t_cnm;
skip_frames = condition_caiman_data.frames_skipped + condition_caiman_data.t_init;
tot_frames = num_frames + numel(skip_frames);
caiman_frames = setdiff([1:tot_frames],skip_frames);

tot_frames = length(trials.stim_type)*opt.trial_length;
for i = 1:num_cells
    lin_idx = zeros(size(temp_coords,1),1);
    
    for t = 1:size(temp_coords,1)
        lin_idx(t) = sub2ind(cnm_dims,temp_coords(t,1),temp_coords(t,2));
    end

    temp_trace = nan(1,tot_frames);
    temp_trace(caiman_frames) =  condition_caiman_data.filt_C(i,1:num_frames);
    temp_trace = fillmissing(temp_trace,'linear');   
    cell_struct(i).filtC = temp_trace;

end
%% check trigger and target cell traces
sens_stim_frames = condition_caiman_data.sensory_stim_frames+condition_caiman_data.t_init;
plot_frame_range = double([condition_caiman_data.t_init sens_stim_frames(tot_num_trials)+opt.trial_length]);
figure('name','online trigger cell traces','units','normalized','outerposition',[0 0 1 1]); hold on
plot_offset = 5;
cell_count = 1;

trigger_idx = baseline_output.trigger_idx;
this_idx = trigger_idx;

xlim(plot_frame_range)
[~,sorted_idx] = sort(online_w(this_idx));
online_tex_auc = cell2mat({cell_struct(:).('correct_stimAUC_zscore')});
[~,sorted_idx] = sort(online_tex_auc(this_idx));

this_idx = this_idx(sorted_idx);
% for i = 1:length(this_idx)
%     ii = this_idx(i);
%     cell_count = cell_count+1;
%     this_trace = condition_cell_struct(ii).filtC;
%     this_trace = zscore(this_trace);
%     plot(this_trace+cell_count*plot_offset,'black');
%     this_cell_color = 'black';
%     if ~isempty(intersect(ii,baseline_output.target_ensembles{1}))
%         this_cell_color = trial_color.stim1;
%     elseif ~isempty(intersect(ii,baseline_output.target_ensembles{2}))
%         this_cell_color = trial_color.stim2;
%     end
%     
%     text(plot_frame_range(1),cell_count*plot_offset,['Cell ' num2str(ii) ', W ' num2str(online_w(ii))], 'horizontalalignment','right', 'color',this_cell_color)
%     
% end
% 
plot_traces = cell2mat({condition_cell_struct(this_idx).filtC}');
plot_traces = plot_traces./max(plot_traces,[],2);
imagesc(plot_traces)
colormap(flipud(gray))
% 
% set(gca,'ytick',[])
% ylim([-plot_offset plot_offset*(cell_count+1)])

% mark trial times
for i = 1:tot_num_trials
    this_stim = trials.stim_type(i);
    this_color = trial_color.(['stim' num2str(this_stim)]);
    p(this_stim)= plot([sens_stim_frames(i) sens_stim_frames(i)]+opt.stimon_frame,ylim,'color',this_color); % trial-on + roughly first touch
%     p(this_stim) = plot([sens_stim_frames(i) sens_stim_frames(i)]+opt.gocue_frame,ylim,'color',this_color); % go-cue
end

photo_stim_frames =  condition_caiman_data.online_photo_frames + condition_caiman_data.t_init;
photo_stim_frames(photo_stim_frames>sens_stim_frames(end)+opt.trial_length)=[];

for i = 1:numel(photo_stim_frames)
    p(3) = plot([photo_stim_frames(i) photo_stim_frames(i)],ylim,'color',trial_color.photostim,'linestyle',':');
end
% 
% for i = 1:numel(dummy_photo_stim_frames)
%     plot([dummy_photo_stim_frames(i) dummy_photo_stim_frames(i)],ylim,'color',trial_color.photostim,'linestyle',':')
% end
% hold on;plot(10.*online_traj./max(abs(online_traj)),'color',[.5 .5 .5])
% text(plot_frame_range(1),0,['Pop. trajectory'], 'horizontalalignment','right', 'color','black')
% xlim(plot_frame_range)
% xlabel('Frames')
% set(gca,'xticklabel',[]);
legend(p,'nogo-tex','go-tex','photo','location','northeastoutside','box','off')
export_fig([fig_save_path filesep 'OnlineConditionTraces' '.pdf'])


%% plot online recorded trajectory
plot_fds = {'stim_1_photo','stim_1_dummyphoto','stim_1_nonphoto'...
    'stim_2_photo','stim_2_dummyphoto','stim_2_nonphoto'}; % from test sta struct
plot_fds = {'stim_1_photo','stim_1_nonphoto'...
    'stim_2_photo','stim_2_nonphoto'}; % from test sta struct
plot_fds = {'stim_1_photo_online','stim_1_nonphoto_online'...
    'stim_2_photo_online','stim_2_nonphoto_online'}; % from test sta struct

plot_num_cols = 4;
plot_num_rows = 2;
ylimit = [-0.5 0.5];
photo_enable_time = [tex_output_struct.pop_params.frames_enable_trigger-opt.gocue_frame]./opt.frame_rate;

figure('name','ensemble traces','position',[100 100 800 600])

for f = 1:numel(plot_fds)
    this_fd = plot_fds{f};
    
    subplot(plot_num_rows,plot_num_cols,f)
    hold on
    this_trials = trial_indices.(this_fd);
    this_traces = traj_struct.(this_fd);
    this_xticks = ([1:opt.trial_length]-opt.sta_gocue_frame)./opt.frame_rate;
%     shadedErrorBar(this_xticks,nanmean(this_traces,1),...
%         nanstd(this_traces,[],1),{'color',[.5 .5 .5],'linewidth',2},0.1);
    shadedErrorBar(this_xticks,median(this_traces,1),[quantile(this_traces,0.75)-median(this_traces,1);...
        median(this_traces,1)-quantile(this_traces,0.25)],{'color',[.5 .5 .5],'linewidth',2},0.5)
    xlabel('Time from go-cue (sec)')
    ylim(ylimit)
    xlim([-4,4])
    title(strrep(this_fd,'_',' '))
    plot([0 0],ylim,':','color',[.5 .5 .5],'linewidth',2) % go-cue
    plot(xlim,[0 0],'black');
    plot(photo_enable_time(1).*[1 1],ylim,':','color',trial_color.photostim,'linewidth',2) % go-cue
    ylabel('Projection')
    axis square
    
    
    subplot(plot_num_rows,plot_num_cols,f+plot_num_cols)
    hold on
    this_trials = trial_indices.(this_fd);
    this_traces = traj_struct.(this_fd);
    this_xticks = ([1:opt.trial_length]-opt.sta_gocue_frame)./opt.frame_rate;
%     shadedErrorBar(this_xticks,mean(this_traces,1),...
%         std(this_traces,[],1),{'color',[.5 .5 .5],'linewidth',2},0.1);
    shadedErrorBar(this_xticks,median(this_traces,1),[quantile(this_traces,0.75)-median(this_traces,1);...
        median(this_traces,1)-quantile(this_traces,0.25)],{'color',[.5 .5 .5],'linewidth',2},0.5)
    plot(repmat(this_xticks,[size(this_traces,1),1])',this_traces','color',[.5 .5 .5])
    
    % mark photostim frames
    if any(isnan(first_photo_frames1(this_trials))~=1)
        this_photo_frames = first_photo_frames1(this_trials)+opt.sta_baseline_frames+1;
        this_photo_times =this_xticks( this_photo_frames );
        scatter(this_photo_times,arrayfun(@(x)this_traces(x,this_photo_frames(x)),1:size(this_traces,1)),'filled','v','markerfacecolor',trial_color.photostim,'markeredgecolor','none')
    end
    xlabel('Time from go-cue (sec)')
    ylim(ylimit)
    xlim([-1,1])
    title(strrep(this_fd,'_',' '))
    plot([0 0],ylim,':','color',[.5 .5 .5],'linewidth',2) % go-cue
    plot(xlim,[0 0],'black');
    plot(photo_enable_time(1).*[1 1],ylim,':','color',trial_color.photostim,'linewidth',2) % go-cue
    ylabel('Projection')
    axis square
    
end

% opt.trial_color = trial_color;
% plot_pop_vectors(traj_struct,plot_fds,1,opt,...
%     'noise_thresh',0,'plot_ylabel','Projection',...
%     'ylimit',ylimit,'IF_MEDIAN',1,...
%     'plot_num_cols',plot_num_cols,'IF_PLOT_RAW_ONLY',0)
suptitle('Closed-loop condition trials, online trajectory')
% 
% control_fds_of_interest = {'stim_5_photostim_1','stim_5_photostim_2','stim_5_nonphotostim'};
% plot_pop_vectors(traj_struct,control_fds_of_interest,1,opt,...
%     'noise_thresh',thresh_sd,'plot_ylabel','Projection',...
%     'plot_num_cols',3,'IF_PLOT_RAW_ONLY',0)
% suptitle('Catch condition trials, online trajectory')
fig=gcf;
set(fig,'PaperPositionMode','auto');         
set(fig,'PaperOrientation','landscape');
% export_fig([fig_save_path filesep 'OnlineTrajectoryPDF' '.pdf'])
export_fig 'D:\TextureData\data\gonogo\\cb291\20200309\\figures\OnlineTrajectoryPDF.pdf'


%% simulation using post-hoc decoder
% get stas


%% target response during photoexcitabitlity test -  didnt manage to align 
photo_opt.N = 1; % auc zscore threshold
photo_opt.sta_amp_thresh = 0.1; % sta amplitude threshold

photo_opt.pre_exp_frames = 0;
photo_opt.sta_pre_frames = 20;
photo_opt.sta_post_frames = 30;
photo_opt.sta_baseline_frames = 20;
photo_opt.bs_frame_range = photo_opt.sta_pre_frames+[(-photo_opt.sta_baseline_frames+1):0];
photo_opt.sta_avg_frames = 20; % frames after stim to average for response amp
photo_opt.sta_thresh = 1;
photo_opt.frame_rate = 30; % Hz

%% load CAIMAN data
[caiman_file,caiman_path] = uigetfile('*.mat','Select seq photostim caiman data');
disp(['Loaded file :',fullfile(caiman_path,caiman_file)])
caiman_data = load(fullfile(caiman_path,caiman_file));
[file,path] = uigetfile('*.mat','Select cell identity data',caiman_path);
disp(['Loaded file :',fullfile(path,file)])
cell_identity_data = load(fullfile(path,file));
target_cell_idx = cell_identity_data.output.target_idx; % this is idx in cell_struct of OnlineProcTexture

%% make data struct
[cnm_struct,cnm_image,num_comp,cnm_dims,tot_frames] = make_cnm_struct(caiman_data);
% spatial components
com_fov = zeros(cnm_dims);
binary_fov = zeros(cnm_dims);
for i = 1:num_comp
    com_fov = com_fov+cnm_struct(i).shape;
end
accepted_idx = caiman_data.accepted_idx+1;
num_cells = numel(accepted_idx);

%% make sta based on photostim sequence
photo_opt.frames_with_photo = double(caiman_data.photoDuration/opt.frame_rate+1); % discard skipped frames with photostim when calculating photo sta amplitude
cell_struct = struct();
photo_stim_frames =  caiman_data.online_photo_frames + caiman_data.t_init;
photo_stim_frames(photo_stim_frames>caiman_data.t_cnm)=[];
temp_photo_sequence_idx = caiman_data.photo_sequence_idx+1; % this is idx in pyRTAOI fixedTargets
temp_photo_sequence_idx = temp_photo_sequence_idx(1:length(photo_stim_frames));
photo_sequence_cell_idx = nan(1,length(photo_stim_frames));
target_cnm_idx = nan(size(target_cell_idx));
for t = 1:numel(target_cell_idx)
    photo_sequence_cell_idx(temp_photo_sequence_idx == t) = target_cell_idx(t); %  this is idx in cell_struct
end

[cell_struct] = make_cell_photostim_sta_struct(cnm_struct,cell_struct,accepted_idx,photo_stim_frames,photo_sequence_cell_idx,photo_opt);

for t = 1:numel(target_cell_idx)
    target_cnm_idx(t) = cell_struct(target_cell_idx(t)).cnm_idx; % this is idx in cnm_struct
end
disp('got photo sta')
%% plot target photoresponse traces
photo_sta_traces{1} = cell2mat( {cell_struct(cell_idx_struct.photo_stim1).sta_trace});
photo_sta_traces{2} = cell2mat( {cell_struct(cell_idx_struct.photo_stim2).sta_trace});
figure
subplot(1,2,1)
plot(photo_sta_traces{1})
subplot(1,2,2)
plot(photo_sta_traces{2})

%% projection on decoder - need to get raw online trace, skipped
decod_proj = struct();
for p = 1:numel(plot_fds)
    fd = plot_fds{p};
    this_traces = arrayfun(@(x)x.(fd).*x.weight,cell_struct,'un',0);
    this_traces = sum(cat(3,this_traces{:}),3)+tex_output_struct.pop_params.thresh;
    decod_proj.(fd) = this_traces;
end

figure
hold on
for f = 1:numel(plot_fds)
    fd = plot_fds{f};
    this_traces = decod_proj.(fd)';
    shadedErrorBar(this_xticks,mean(this_traces,1), std(this_traces,[],1),{'color',trial_color.(fd),'linewidth',2},0.1);
end
ylabel(['Decoder Projection'])
xlabel('Time from go-cue (sec)')
xlim([-5,4])
plot([0 0],ylim,':','color',[.5 .5 .5],'linewidth',2) % go-cue
window_enable_photostim = [[tex_output_struct.pop_params.frames_enable_trigger]+opt.sta_baseline_frames- opt.sta_gocue_frame]./opt.frame_rate;
plot(window_enable_photostim,[0 0],'color',trial_color.photostim,'linewidth',2) % sta average window (withold, 1sec before gocue)

%% === control session ==
%% STA of online trajectory
online_traj = ctr_caiman_data.online_traj;
online_traj = online_traj./max(abs(online_traj));
num_frames = ctr_caiman_data.t_cnm;
skip_frames = ctr_caiman_data.frames_skipped + ctr_caiman_data.t_init+1;
tot_frames = num_frames + numel(skip_frames);
caiman_frames = setdiff([1:tot_frames],skip_frames);

traj_struct = struct();
temp_trace = nan(1,tot_frames);
temp_trace(caiman_frames) = online_traj(1,1:num_frames);
%  online_traj = fillmissing(temp_trace,'linear');
online_traj = temp_trace;
sens_stim_frames = ctr_caiman_data.sensory_stim_frames+ctr_caiman_data.t_init;
this_sens_stim_frames = sens_stim_frames+opt.gocue_frame;
[~,~,~,traj_struct.sta_traces,~,~,traj_struct.sta_trace] =...
    make_sta_from_traces(online_traj,this_sens_stim_frames,opt.sta_pre_frames,opt.sta_post_frames,1:opt.sta_baseline_frames);


