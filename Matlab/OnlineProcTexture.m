% generate sta from pyRTAOI result file after a sensory stim movie
% note: rois idx doesn't match with what's shown in pyRTAOI image window
% run this to quickly find sensory-opsin-positve cells  

% ZZ 2018

% adapted from OnlineProcVisual

% TO DO:
% concatinate multiple sessions
% check different epochs; 
% classify cells baseed on epoch time-average
% classify cells baseed on t-SNE with 'baits'
% check trial average
% condition on catch trials


%% add path - change this for rig
clear all
close all
clc

% ZZ PC
% addpath(genpath('C:\Users\User\Desktop\pyRTAOI-rig\Matlab'));
% cd('C:\Users\User\Desktop\pyRTAOI-rig\Matlab');

% BRUKER1
matlab_set_paths_zz

%% parameters - CHANGE THIS
IF_GO_NOGO = true;
opt.go_stim_type = 2;% miss trials will be taken as incorrect for go stim type
opt.nogo_stim_type = 1;

opt.N = 1.5; % threshold for significant auc
opt.sta_pre_frames = 150; 
opt.sta_post_frames = 110;
opt.sta_baseline_frames = 30; % relative to beginning of sta traces
opt.trial_length = 1+opt.sta_pre_frames+opt.sta_post_frames;

opt.sta_stim_frame = 90; % roughly when first contact occurs
opt.gocue_frame = 120; % relative to trial start
opt.stimon_frame = 90; % relative to trial start
opt.end_rw_frame = 180; % end of response window

% time window for first lick, 
% will not be taken as miss trial if first lick occur during withold window
opt.rw_win_end_sec = 5;
opt.withold_win_start_sec = 3;
opt.rw_start_sec = 3.6; %trials with first lick before this will not be used for analysis; 4 is beginning of response window in pybehav 

% frame indices relative to sta trace
opt.sta_gocue_frame = opt.sta_pre_frames;
opt.sta_trialon_frame = opt.sta_pre_frames-opt.gocue_frame;
opt.sta_avg_frames = [-30:1:-5]+opt.sta_gocue_frame; % 1 sec before go-cue
opt.sta_peak_search_range =  [-60:1:0]+opt.sta_gocue_frame;


opt.withold_frames_adj = [60 120]; 

opt.offcell_avg_frames = 30;
opt.sta_amp_thresh = 1;
opt.frame_rate = 30;

opt.flag_use_peak = false; % if false then will use average to get roc

opt.correct_trial_only = true; % only use correct trials to get tunning

% select which online trace to use
opt.IF_USE_FILTC = true;

% select which dim reduction to use
opt.method = 'dpca'; % fa, dpca, cd

% discard trials after this
opt.discard_trials_after = [];
opt.discard_trials_before = [20];

% define epochs (relative to trial start )
opt.epoch_names = {'pre','early','late'};
opt.epoch_frames = {[1:45],[60:90],[90:120]};

% init color
[trial_color] = online_tex_init_color();
%% load CAIMAN data
[caiman_file,caiman_path] = uigetfile('*.mat','Select texture caiman data');
disp(['Loaded file :',fullfile(caiman_path,caiman_file)])
caiman_data = load(fullfile(caiman_path,caiman_file)); 

%% load CNN results
FLAG_CNN_RESULT_LOADED = false;
try
[cnn_file,cnn_path] = uigetfile([caiman_path '*.mat'],'Select CNN result');
disp(['Loaded file :',fullfile(cnn_path,cnn_file)])
cnn_data = load(fullfile(cnn_path,cnn_file)); 
cnn_thresh = cnn_data.CNN_thresh;
cnn_predictions = double(cnn_data.CNN_predictions);
FLAG_CNN_RESULT_LOADED = true;
catch
    disp('not using cnn results')
    cnn_thresh = 0;
    cnn_predictions = ones(1,caiman_data.cnm_N);
end
%% load Pybehavior data
try
[pb_file,pb_path] = uigetfile('*.mat','Select pybehavior data');
disp(['Loaded file :',fullfile(pb_path,pb_file)])
behavior_data =  load(fullfile(pb_path,pb_file)) ;
trials = make_trials_struct(behavior_data);
FLAG_PYBEHAV_LOADED = true;
catch
    FLAG_PYBEHAV_LOADED = false;
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
opt.exp_name = strrep(caiman_file,'.mat','');

%% setup save path and name
opt.output_path = caiman_path;
opt.exp_name = strrep(caiman_file,'.mat','proc');

%% check if trial order matches in behavior and caiman file
if FLAG_PYBEHAV_LOADED
    if isequal(trials.stim_type,caiman_data.trialOrder)
        disp('data matched')
    else
        warning('trial order mismatch!')
    end
else
    warning('no pybehav data loaded')
end
%% == ORGANISE TRIAL TYPES ==
tex_stim_frames = {};
if(~isempty(caiman_data.stim_frames_caiman))
    sens_stim_frames = caiman_data.sensory_stim_frames+caiman_data.t_init;
else
    sens_stim_frames = [];
end

% discard sens stim frames beyond tottime (deal with senarios when movie is aborted early)

discard_trial_idx = find(sens_stim_frames>caiman_data.t_cnm-opt.sta_post_frames-opt.sta_pre_frames);

if ~isempty(opt.discard_trials_after)
    discard_trial_idx = [discard_trial_idx, opt.discard_trials_after:numel(sens_stim_frames)];
end
if ~isempty(opt.discard_trials_before)
    discard_trial_idx = [discard_trial_idx, 1:opt.discard_trials_before];
end
discard_trial_idx = unique(discard_trial_idx);
sens_stim_frames(discard_trial_idx) = [];
if FLAG_PYBEHAV_LOADED
    trials.miss(trials.firstlick<opt.rw_start_sec) = 0; % discard early response trials from miss trials
    trials.correct(trials.firstlick<opt.rw_start_sec) = 0; % discard early response trials 
    trials.incorrect(trials.firstlick<opt.rw_start_sec) = 0; % discard early response trials
    trials.fa(trials.firstlick<opt.rw_start_sec) = 0; % discard early response trials

    fds = fields(trials);
    try
        for f =1:numel(fds)
            trials.(fds{f})(discard_trial_idx) = [];
        end
    end
end
% trial type
trialOrder = caiman_data.trialOrder;
try
trialOrder(discard_trial_idx) = [];
catch
    disp('no discarded trials')
end
trialTypes = unique(trialOrder);
num_stim_type = length(unique(trialOrder)); % orientations/texture

if ~ FLAG_PYBEHAV_LOADED
    % make a dummy trial struct (setting all trials to correct)
    trials.correct = ones(size(sens_stim_frames));
    trials.incorrect = zeros(size(sens_stim_frames));
    trials.stim_type = trialOrder;
    trials.cheated = zeros(size(sens_stim_frames));
    trials.miss = zeros(size(sens_stim_frames));
end

for i = 1:numel(trialTypes)
    tex_stim_frames{i} = sens_stim_frames(trialOrder == trialTypes(i));
end

% roi indices
opsin_positive = caiman_data.opsin_positive;
accepted_idx = caiman_data.accepted_idx+1;
opsin_positive_idx = accepted_idx(opsin_positive>0);

% color trial by stim type
% hsv_lut = colormap(hsv);
% hsv_lut = hsv_lut(2:end-3,:);
% close
% indices = round(linspace(1,size(hsv_lut,1),num_stim_type));

% opt.trial_color = zeros(numel(trialOrder),3);
% for t = 1:numel(trialOrder)
%     opt.trial_color(t,:) = trial_color.(['stim' num2str(trialOrder(t))]);
% end

num_trials = numel(sens_stim_frames);
opt.type_color = [trial_color.('stim1');trial_color.('stim2')];
opt.trial_color = trial_color;
%% == ORGANISE IMAGING DATA ==
%% make cnm data structure
cnm_struct = struct();
cnm_dims = double(caiman_data.cnm_dims);
cnm_image = reshape(caiman_data.cnm_b,cnm_dims);
cnm_A = full(caiman_data.cnm_A);
num_frames = caiman_data.t_cnm;

num_comp = size(cnm_A,2);
comp_shape =[cnm_dims(1),cnm_dims(2),num_comp];
cm = com(sparse(double(cnm_A)),cnm_dims(1),cnm_dims(2)); % center of mass


if ~isempty(caiman_data.frames_skipped)
    skip_frames = caiman_data.frames_skipped + caiman_data.t_init;
    tot_frames = num_frames + numel(skip_frames);
    caiman_frames = setdiff([1:tot_frames],skip_frames);
else
    caiman_frames = 1:num_frames;
    tot_frames = num_frames;
end

for i = 1:num_comp
    cnm_struct(i).shape = reshape(cnm_A(:,i),cnm_dims);
    cnm_struct(i).onlineC = caiman_data.online_C(i+1,1:num_frames);
    cnm_struct(i).deconvC = caiman_data.cnm_C(i,1:num_frames);
    cnm_struct(i).centroid = cm(i,:);
    
    % set skipped frames to nan
    temp_trace = nan(1,tot_frames);
    temp_trace(caiman_frames) =  caiman_data.cnm_C(i,1:num_frames);
    temp_trace = fillmissing(temp_trace,'linear');
    cnm_struct(i).deconvC_full = temp_trace;
    
    temp_trace = nan(1,tot_frames);
    temp_trace(caiman_frames) =  caiman_data.noisyC(i+1,1:num_frames);
    temp_trace = fillmissing(temp_trace,'linear');
    cnm_struct(i).noisyC_full = temp_trace;
    
    % use go-cue frame as 'stim frame' for sta traces
    cnm_struct(i).stim_frames = sens_stim_frames+opt.gocue_frame;

end

temp_trace = nan(1,tot_frames);
temp_trace(caiman_frames) =  caiman_data.noisyC(1,1:num_frames);
temp_trace = fillmissing(temp_trace,'linear');
backgroundC = temp_trace;

disp('made cnm struct')

% only get accepted cells
accepted_idx = caiman_data.accepted_idx+1;
num_cells = numel(accepted_idx);
%% plot spatial components and save to cell struct
com_fov = zeros(cnm_dims);
accepted_com_fov = zeros(cnm_dims);
for i = 1:num_comp
    com_fov = com_fov+cnm_struct(i).shape;
end

for i = accepted_idx
    accepted_com_fov = com_fov+cnm_struct(i).shape;
end


cnm_plot_options = CNMFSetParms;
cnm_plot_options.roi_color = [colormap(lines);colormap(lines);colormap(lines)];
close

figure('name','fov','position',[100 100 1200 500])
subplot(1,3,1)
imagesc(com_fov)
colormap(gray)
axis square
title('Detected ROIs')

subplot(1,3,2)
imagesc(accepted_com_fov)
colormap(gray)
axis square
title('Accepted ROIs')


subplot(1,3,3)
[CC,jsf] = plot_contours(sparse(double(cnm_A)),cnm_image,cnm_plot_options,1,[],[],[1 1 1]);
colormap(gray)
axis square
title('GCaMP')

cell_struct = struct();
for i = 1:num_cells
    this_idx = accepted_idx(i);
    temp_coords = jsf(this_idx).coordinates;
    lin_idx = zeros(size(temp_coords,1),1);
   
    for t = 1:size(temp_coords,1)
        lin_idx(t) = sub2ind(cnm_dims,temp_coords(t,1),temp_coords(t,2));
    end
    cell_struct(i).contour = CC{this_idx};
    cell_struct(i).lin_coords = lin_idx;
    cell_struct(i).coordinates = jsf(this_idx).coordinates;
    cell_struct(i).pix_values = jsf(this_idx).values;
    cell_struct(i).centroid = jsf(this_idx).centroid;
    cell_struct(i).opsin_positive = 0;
    cell_struct(i).cnm_idx = this_idx;
    cell_struct(i).jsf = jsf(this_idx);
    cell_struct(i).cnn_prediction = cnn_predictions(this_idx);
    cell_struct(i).filtC = caiman_data.filt_C(i,1:num_frames);

    if(~isempty(find(opsin_positive_idx==i)))
         cell_struct(i).opsin_positive = 1;
    end
end

%% plot traces
figure('name','traces','position',[100 100 1600 800]); hold on
plot_offset = 5;
stim_cell_count = 1;
non_stim_cell_count = 1;

for i = 1:num_cells
    this_cell_trace = zscore(cnm_struct(cell_struct(i).cnm_idx).onlineC);
    plot(this_cell_trace+i*plot_offset,'color','black','linewidth',1.5)
    stim_cell_count = stim_cell_count+1;
end
xlabel('Frames')
ylabel('ROI index')
yticks([ 1:num_comp].*plot_offset)
yticklabels(1:num_comp)
xlim([caiman_data.t_init tot_frames])
ylim([0 num_comp].*plot_offset+5)

% background
plot(backgroundC,'color',[.5 .5 .5],'linestyle',':')

for i = 1:numel(sens_stim_frames)
    this_color = trial_color.(['stim' num2str(trials.stim_type(i) )]);
    plot([sens_stim_frames(i) sens_stim_frames(i)],ylim,'color',this_color,'linestyle',':') % trial-on 
    plot([sens_stim_frames(i) sens_stim_frames(i)]+opt.gocue_frame,ylim,'color',this_color) % go-cue
end
export_fig([fig_save_path filesep 'OnlineTraces_' strrep(caiman_file,'.mat','.png')])

%% Get stim triggered average 
for i = 1:num_cells
    this_cell_trace = zscore( cnm_struct(cell_struct(i).cnm_idx).deconvC_full);
    
    if opt.IF_USE_FILTC
        this_online_trace = cell_struct(i).filtC; % online trace med filtered
    else
        this_online_trace = cnm_struct(cell_struct(i).cnm_idx).onlineC_full; % online trace without median filter
    end

    
    this_num_trials = numel(cnm_struct(cell_struct(i).cnm_idx).stim_frames );
    this_sens_stim_frames =  cnm_struct(cell_struct(i).cnm_idx).stim_frames;
    cell_struct(i).num_trials = this_num_trials;
    cell_struct(i).is_sensory = 0;
    cell_struct(i).is_offcell = 0;
    cell_struct(i).pref_orient = [];
    cell_struct(i).sta_amp = 0;
    cell_struct(i).sta_traces = [];
    cell_struct(i).sta_trace = [];
    cell_struct(i).accepted = 0;
    cell_struct(i).pref_orient = nan;
    if(~isempty(find(accepted_idx==i)))
        cell_struct(i).accepted = 1;
    end
    
    if(this_num_trials>0)
        % average across all stim types
        % using df 
%         [~,~,~,~,~,cell_struct(i).sta_traces,cell_struct(i).sta_trace] = make_sta_from_traces(this_cell_trace,this_sens_stim_frames,opt.sta_pre_frames,opt.sta_post_frames,1:opt.sta_baseline_frames);
          % using raw f
        [~,~,~,cell_struct(i).sta_traces,~,~,cell_struct(i).sta_trace] = make_sta_from_traces(this_cell_trace,this_sens_stim_frames,opt.sta_pre_frames,opt.sta_post_frames,1:opt.sta_baseline_frames);
        [~,~,~,cell_struct(i).raw_sta_traces,~,~,cell_struct(i).raw_sta_trace] = make_sta_from_traces(this_online_trace,this_sens_stim_frames,opt.sta_pre_frames,opt.sta_post_frames,1:opt.sta_baseline_frames);

            cell_struct(i).sta_amp = mean(cell_struct(i).sta_trace(opt.sta_avg_frames));

        if  cell_struct(i).sta_amp > opt.sta_amp_thresh
            cell_struct(i).is_sensory = 1;
        end        
    end
    
end

% sta of backgound component (for alignment check)
[~,~,~,~,~,bg_sta_traces,bg_sta_trace] =...
    make_sta_from_traces(backgroundC,sens_stim_frames+opt.sta_gocue_frame ,opt.sta_pre_frames,opt.sta_post_frames,1:opt.sta_baseline_frames);
disp('got cell_struct sta')
%% Normalise traces to baseline and std
% get baseline and std from the first sets of easy trials
% this is implemented for online - shouldnt need to do again here?
easy_trial_idx = 1:10;
cell_bs = cell2mat(arrayfun(@(x)mean(reshape(cell_struct(x).raw_sta_traces(easy_trial_idx,1:opt.sta_baseline_frames),[1,numel(easy_trial_idx)*opt.sta_baseline_frames]),2),1:num_cells,'UniformOutput', false));
cell_sd = cell2mat(arrayfun(@(x)std(reshape(cell_struct(x).raw_sta_traces(easy_trial_idx,:),[1,numel(easy_trial_idx)*opt.trial_length])),1:num_cells,'UniformOutput', false));
for i = 1:num_cells   
    cell_struct(i).raw_sta_traces = (cell_struct(i).raw_sta_traces - cell_bs(i));
    cell_struct(i).raw_sta_trace = (cell_struct(i).raw_sta_trace - cell_bs(i));
end
disp('normalised cell_struct raw_sta_traces')


%% Get trial types
trial_indices = struct(); % % get trialtype-outcome indices
for i = 1:num_stim_type
    trial_indices.(['stim_' num2str(i) '_correct']) = find(trials.correct==1&trials.stim_type==i&trials.cheated==0);
    trial_indices.(['stim_' num2str(i) '_incorrect']) = find(trials.incorrect==1&trials.stim_type==i&trials.cheated==0);
    trial_indices.(['stim_' num2str(i) '_correct'])(trial_indices.(['stim_' num2str(i) '_correct'])>num_trials)=[];
    trial_indices.(['stim_' num2str(i) '_incorrect'])(trial_indices.(['stim_' num2str(i) '_incorrect'])>num_trials)=[];
    
end
if IF_GO_NOGO %overide incorrect trials by miss trials for go-stim
    trial_indices.(['stim_' num2str(opt.go_stim_type) '_incorrect']) = find(trials.miss==1&trials.stim_type==opt.go_stim_type&trials.cheated==0& ...
        (trials.firstlick>opt.rw_win_end_sec|trials.firstlick<opt.withold_win_start_sec|isnan(trials.firstlick)));
    % no-go stim type correct trials should exclude trials with licks in
    % withold window (use fa as incorrect)
    trial_indices.(['stim_' num2str(opt.nogo_stim_type) '_correct']) = find(trials.correct==1&trials.fa==0&trials.stim_type==opt.nogo_stim_type&trials.cheated==0);
    trial_indices.(['stim_' num2str(opt.nogo_stim_type) '_incorrect']) = find(trials.fa==1&trials.stim_type==opt.nogo_stim_type&trials.cheated==0);

end
try
% catch trials
    trial_indices.(['stim_5_miss']) = find((trials.stim_type==3|trials.stim_type==4)&trials.miss==1); 
    trial_indices.(['stim_5_lick']) = find((trials.stim_type==3|trials.stim_type==4)&trials.miss==0); 

end

[dp,hr,fa] = get_gonogo_baseline_performance(trial_indices);


%% Sort STAs for different trial types
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


%% == GET CELL IDENTITY ==
% stimulus AUC calculated from correct trials
% takes long, consider make it simple 
peak_frame_range = opt.sta_peak_search_range;
avg_frame_range = opt.sta_avg_frames;
[cell_struct] = get_cell_auc(cell_struct,{'stim_1_correct','stim_2_correct'},'correct_stimAUC',opt);
[cell_struct] = get_cell_active_auc(cell_struct,{'stim_1_correct','stim_2_correct'},opt);

correct_stimulusAUC_zscore = extractfield(cell_struct,'correct_stimAUC_zscore');
stim_1_correct_offcell = extractfield(cell_struct,'stim_1_correct_offcell');
stim_2_correct_offcell = extractfield(cell_struct,'stim_2_correct_offcell');


correct_texAUC = extractfield(cell_struct,'correct_stimAUC');
disp('got correct trials auc')

% try
%     [cell_struct] = get_cell_auc(cell_struct,{'stim_1_incorrect','stim_2_incorrect'},'incorrect_stimAUC',opt);
%     FLAG_GOT_CHOICE_AUC = true;
% catch
%     disp('getting incorrect trials error')
     FLAG_GOT_CHOICE_AUC = false;
% end

disp('got stim auc')
%% get cell_idx_struct
cell_idx_struct = struct();
cell_idx_struct.all = 1:num_cells;
cell_idx_struct.all_tex2 = find(correct_stimulusAUC_zscore>opt.N); % cells prefering texture1 in correct trials
cell_idx_struct.all_tex1 = find(correct_stimulusAUC_zscore<-opt.N); % cells prefering texture2 in correct trials

% exclude cells with negative sensory response from target
cell_idx_struct.tex2 = find(correct_stimulusAUC_zscore>opt.N&stim_1_correct_offcell==0); % cells prefering texture2 in correct trials
cell_idx_struct.tex1 = find(correct_stimulusAUC_zscore<-opt.N&stim_2_correct_offcell==0); % cells prefering texture2 in correct trials

cell_idx_struct.tex = unique([find(correct_stimulusAUC_zscore>opt.N), find(correct_stimulusAUC_zscore<-opt.N)]);

if FLAG_GOT_CHOICE_AUC
incorrect_stimulusAUC_zscore = extractfield(cell_struct,'incorrect_stimAUC_zscore');
incorrect_texAUC = extractfield(cell_struct,'incorrect_stimAUC');
cell_idx_struct.port2 = find(correct_stimulusAUC_zscore>opt.N & incorrect_stimulusAUC_zscore<-opt.N); % cells prefering texture1 in correct trials
cell_idx_struct.port1 = find(correct_stimulusAUC_zscore<-opt.N& incorrect_stimulusAUC_zscore>opt.N); % cells prefering texture2 in correct trials
cell_idx_struct.port = unique([cell_idx_struct.port1,cell_idx_struct.port2]);
cell_idx_struct.relevant = unique([cell_idx_struct.port,cell_idx_struct.tex]);
else
    warning('no choice auc')
end
for c = 1:num_cells
    cell_struct(c).is_tuned = 0;
    cell_struct(c).pref = 0;
    cell_struct(c).is_offcell = 0;
    if ismember(c,cell_idx_struct.tex)
        cell_struct(c).is_tuned = 1;
        if ismember(c,cell_idx_struct.all_tex1)
            cell_struct(c).pref = 1;
             cell_struct(c).is_offcell =  cell_struct(c).stim_1_correct_offcell;
        else
            cell_struct(c).pref = 2;
            cell_struct(c).is_offcell =  cell_struct(c).stim_2_correct_offcell;
        end
    end
end
disp('got cell identity')
%% check cells with low snr or bad shape (cnn prediction score)
cell_snr = cell_bs./cell_sd.^2;
cell_snr(cell_sd<0.001) = 0;
cell_size = cell2mat(arrayfun(@(x)numel(x.pix_values),cell_struct,'UniformOutput', false));
for  c = 1:num_cells
    cell_struct(c).snr = cell_snr(c);
    cell_struct(c).size = cell_size(c);
end


%% Discard cells of bad shapes from trigger and targets pool (optional)
% not using cells with low cnn score as trigger or targets
% will not do anything if no cnn result is loaded
% cnn_thresh = 0.05; % change this to higher value to exlude more dendrites
cell_idx_struct.cnn_above_thresh = find(cnn_predictions(accepted_idx)>cnn_thresh);
cell_idx_struct = structfun(@(x)intersect(x,cell_idx_struct.cnn_above_thresh),cell_idx_struct,'UniformOutput',false);
disp(['discarded cells with cnn prediction soore <' num2str(cnn_thresh)])
%% select cell identity for readout and stimulation
% check ROI quality - increase cnn_thresh and update cell_idx_struct
opt.target_idx_fd = {'all_tex1','all_tex2'};
opt.trigger_idx_fd = 'tex';
opt.fov_size = double(cnm_dims);
opt.ds_factor = caiman_data.ds_factor;

% check trigger and target cells on FOV
figure('name','trigger targets check','position',[100 100 1600 800]);
ax = subplot(1,3,1);
plot_value_in_rois( cell_struct, 'correct_stimAUC_zscore',[256 256],ax,'IF_NORM_PIX',0,'IF_CONTOUR',0);
set(gca,'Ydir','reverse')
title('Texture selectivity (auc zscore)')

% mark texture cells on FOV
mark_cells_on_fov(cell_struct,cell_idx_struct.(opt.trigger_idx_fd),[.5 .5 .5],'MarkerSize',300,'Linewidth',2)
for f = 1:numel(opt.target_idx_fd)
    fd = opt.target_idx_fd{f};
    mark_cells_on_fov(cell_struct,cell_idx_struct.(fd),trial_color.(['stim' num2str(f)]),'MarkerSize',500,'Linewidth',2)
end

ax = subplot(1,3,2);
plot_value_in_rois( cell_struct, 'cnn_prediction',[256 256],ax,'IF_NORM_PIX',0,'IF_CONTOUR',0);
set(gca,'Ydir','reverse')
title('CNN prediction score')

% mark texture cells on FOV
mark_cells_on_fov(cell_struct,cell_idx_struct.(opt.trigger_idx_fd),[.5 .5 .5],'MarkerSize',300,'Linewidth',2)
for f = 1:numel(opt.target_idx_fd)
    fd = opt.target_idx_fd{f};
    mark_cells_on_fov(cell_struct,cell_idx_struct.(fd),trial_color.(['stim' num2str(f)]),'MarkerSize',500,'Linewidth',2)
end

ax = subplot(1,3,3);

plot_value_in_rois( cell_struct, 'snr',[256 256],ax,'IF_NORM_PIX',0,'IF_CONTOUR',0);
set(gca,'Ydir','reverse')
title('Cell SNR')

% mark texture cells on FOV
mark_cells_on_fov(cell_struct,cell_idx_struct.(opt.trigger_idx_fd),[.5 .5 .5],'MarkerSize',300,'Linewidth',2)
for f = 1:numel(opt.target_idx_fd)
    fd = opt.target_idx_fd{f};
    mark_cells_on_fov(cell_struct,cell_idx_struct.(fd),trial_color.(['stim' num2str(f)]),'MarkerSize',500,'Linewidth',2)
end

% export_fig([fig_save_path filesep 'SelCNNonFOV_' strrep(caiman_file,'.mat','.png')])


%% MAKE OUTPUT FILE FOR PYRTAOI PHOTOEXCITABILITY TEST
[output1,save_time1] = generate_cell_idx_file(cell_struct,cell_idx_struct,[],opt);
disp(['No. photo required: ' num2str(length(output1.target_idx)*15)])
%% Plot STAs traces for all components (save figures) - COULD BE REALLY SLOW
figure('name','condition sta traces (raw)','units','normalized','outerposition',[0 0 1 1])
num_plot_cols = 10;
num_plot_rows = ceil(num_cells/num_plot_cols);
for i = 1:num_cells
    subtightplot(num_plot_rows,num_plot_cols,i)
    hold on
    % correct trials
    for j = 1:num_stim_type
        plot(cell_struct(i).(['stim_' num2str(j) '_correct']),'color',trial_color.(['correct_stim' num2str(j)]),'linewidth',1)
        plot(cell_struct(i).(['stim_' num2str(j) '_incorrect']),'color',trial_color.(['incorrect_stim' num2str(j)]),'linewidth',1)
    end
    
    xlim([0 opt.sta_pre_frames+opt.sta_post_frames+1])
    axis square
    
    text(1,1,['Cell ' num2str(i) ' (ROI ' num2str(cell_struct(i).cnm_idx) ')'],'units','normalized','color','black','Horizontalalignment','right','VerticalAlignment','top')

    if(~isempty(find(opsin_positive_idx==i)))
        text(1,1,['Cell ' num2str(i)],'units','normalized','color','r','Horizontalalignment','right','VerticalAlignment','top')
    end
    box off
    
    % mark tuned cells
    if( cell_struct(i).is_tuned)
        box on
        this_color = trial_color.(['correct_stim' num2str(cell_struct(i).pref)]);
        set(gca,'XColor',this_color,'YColor',this_color,'linewidth',3)
    end
    
    % mark time
    plot([1,1].*opt.sta_gocue_frame, ylim,'color','black','linestyle',':')
    if ~ opt.flag_use_peak
        plot([avg_frame_range(1),avg_frame_range(end)],[0,0],'color',[.5 .5 .5],'linewidth',2)
    else
        plot([peak_frame_range(1),peak_frame_range(end)],[0,0],'color',[.5 .5 .5],'linewidth',2)
    end
    
    % show auc
    text(0.05,.8,['trial auc ' num2str(cell_struct(i).correct_stimAUC,'%10.2f')],'units','normalized', 'horizontalalignment','left', 'color','black')
    text(0.05,.7,['zscore auc '  num2str(cell_struct(i).correct_stimAUC_zscore,'%10.2f') ],'units','normalized', 'horizontalalignment','left', 'color','black')
    
    % modify x ticks
    xaxisvalues = [0:30:opt.sta_pre_frames+opt.sta_post_frames];
    xticks(xaxisvalues)
    xticklabels(arrayfun(@(x){num2str(x)},(xaxisvalues-opt.sta_trialon_frame)./opt.frame_rate))
    
end
export_fig([fig_save_path filesep 'RawSTATrace_' strrep(caiman_file,'.mat','.png')])
%
figure('name','condition sta traces (errobar)','units','normalized','outerposition',[0 0 1 1])
x_ticks =[0:1:opt.trial_length-1];

for i = 1:num_cells
    subtightplot(num_plot_rows,num_plot_cols,i)
    hold on
    for j = 1:num_stim_type
        this_traces = cell_struct(i).(['stim_' num2str(j) '_correct'])';
        shadedErrorBar(x_ticks,mean(this_traces,1),...
            std(this_traces,[],1),{'color',trial_color.(['correct_stim' num2str(j)]),'linewidth',2},0.1);
        try
        this_traces = cell_struct(i).(['stim_' num2str(j) '_incorrect'])';
        shadedErrorBar(x_ticks,mean(this_traces,1),...
            std(this_traces,[],1),{'color',trial_color.(['incorrect_stim' num2str(j)]),'linewidth',2},0.1);
        end
    end
    
    xlim([0 opt.sta_pre_frames+opt.sta_post_frames+1])
    axis square
    
    text(1,1,['Cell ' num2str(i) ' (ROI ' num2str(cell_struct(i).cnm_idx) ')'],'units','normalized','color','black','Horizontalalignment','right','VerticalAlignment','top')

    if(~isempty(find(opsin_positive_idx==i)))
        text(1,1,['Cell ' num2str(i)],'units','normalized','color','r','Horizontalalignment','right','VerticalAlignment','top')
    end
    box off
    
    % mark tuned cells
    if( cell_struct(i).is_tuned)
        box on
        this_color = trial_color.(['correct_stim' num2str(cell_struct(i).pref)]);
        set(gca,'XColor',this_color,'YColor',this_color,'linewidth',3)
    end
    
    % mark time
    plot([1,1].*opt.sta_gocue_frame, ylim,'color','black','linestyle',':','linewidth',2)
    if ~ opt.flag_use_peak
        plot([avg_frame_range(1),avg_frame_range(end)],[0,0],'color',[.5 .5 .5],'linewidth',2)
    else
        plot([peak_frame_range(1),peak_frame_range(end)],[0,0],'color',[.5 .5 .5],'linewidth',2)
    end
    
    % show auc
    text(0.05,.8,['trial auc ' num2str(cell_struct(i).correct_stimAUC,'%10.2f')],'units','normalized', 'horizontalalignment','left', 'color','black')
    text(0.05,.7,['zscore auc '  num2str(cell_struct(i).correct_stimAUC_zscore,'%10.2f') ],'units','normalized', 'horizontalalignment','left', 'color','black')
    
    % modify x ticks
    xaxisvalues = [0:30:opt.sta_pre_frames+opt.sta_post_frames];
    xticks(xaxisvalues)
    xticklabels(arrayfun(@(x){num2str(x)},(xaxisvalues-opt.sta_trialon_frame)./opt.frame_rate))
    
end
export_fig([fig_save_path filesep 'ShadedSTATrace_' strrep(caiman_file,'.mat','.png')])

%% Plot STAs traces for certain cell types
plot_cell_type = 'tex';
plot_cell_idx = cell_idx_struct.(plot_cell_type);
num_plot_cols = 8;
num_plot_rows = ceil(numel(plot_cell_idx)/num_plot_cols);

plot_cell_idx = sort(plot_cell_idx);
figure('name','tex cell sta traces','units','normalized','outerposition',[0 0 1 1])
for ii = 1:numel(plot_cell_idx)
    i = plot_cell_idx(ii);
    subtightplot(num_plot_rows,num_plot_cols,ii)
    hold on
    % correct trials
    for j = 1:num_stim_type
        plot(raw_cell_struct(i).(['stim_' num2str(j) '_correct']),'color',trial_color.(['correct_stim' num2str(j)]),'linewidth',1)
        plot(raw_cell_struct(i).(['stim_' num2str(j) '_incorrect']),'color',trial_color.(['incorrect_stim' num2str(j)]),'linewidth',1)
    end
    
    xlim([0 opt.sta_pre_frames+opt.sta_post_frames+1])
    axis square
    
    text(1,1,['Cell ' num2str(i) ' (ROI ' num2str(cell_struct(i).cnm_idx) ')'],'units','normalized','color','black','Horizontalalignment','right','VerticalAlignment','top')

    if(~isempty(find(opsin_positive_idx==i)))
        text(1,1,['Cell ' num2str(i)],'units','normalized','color','r','Horizontalalignment','right','VerticalAlignment','top')
    end
    box off
    
    % mark tuned cells
    if( cell_struct(i).is_tuned)
        box on
        this_color = trial_color.(['correct_stim' num2str(cell_struct(i).pref)]);
        set(gca,'XColor',this_color,'YColor',this_color,'linewidth',3)
    end
    
    if( cell_struct(i).is_offcell)
        box on
        this_color = [.5 .5 .5];
        set(gca,'XColor',this_color,'YColor',this_color,'linewidth',1)
    end
    
    % mark time
    plot([1,1].*opt.sta_gocue_frame, ylim,'color','black','linestyle',':')
    if ~ opt.flag_use_peak
        plot([avg_frame_range(1),avg_frame_range(end)],[0,0],'color',[.5 .5 .5],'linewidth',2)
    else
        plot([peak_frame_range(1),peak_frame_range(end)],[0,0],'color',[.5 .5 .5],'linewidth',2)
    end
    
    % show auc
    text(0.05,.8,['trial auc ' num2str(cell_struct(i).correct_stimAUC,'%10.2f')],'units','normalized', 'horizontalalignment','left', 'color','black')
    text(0.05,.7,['zscore auc '  num2str(cell_struct(i).correct_stimAUC_zscore,'%10.2f') ],'units','normalized', 'horizontalalignment','left', 'color','black')
    
    % modify x ticks
    xaxisvalues = [0:30:opt.sta_pre_frames+opt.sta_post_frames];
    xticks(xaxisvalues)
    xticklabels(arrayfun(@(x){num2str(x)},(xaxisvalues-opt.sta_trialon_frame)./opt.frame_rate))
    
%     ylim([-5 20])
end
export_fig([fig_save_path filesep plot_cell_type 'STATrace_' strrep(caiman_file,'.mat','.png')])


%%
figure('name','shaded tex cell sta traces','units','normalized','outerposition',[0 0 1 1])
x_ticks =[0:1:opt.trial_length-1];

for  ii = 1:numel(plot_cell_idx)
    i = plot_cell_idx(ii);
    subtightplot(num_plot_rows,num_plot_cols,ii)
    hold on
    for j = 1:num_stim_type
        try
        this_traces = raw_cell_struct(i).(['stim_' num2str(j) '_correct'])';
        shadedErrorBar(x_ticks,mean(this_traces,1),...
            std(this_traces,[],1),{'color',trial_color.(['correct_stim' num2str(j)]),'linewidth',2},0.1);
        end
        try
        this_traces = raw_cell_struct(i).(['stim_' num2str(j) '_incorrect'])';
        shadedErrorBar(x_ticks,mean(this_traces,1),...
            std(this_traces,[],1),{'color',trial_color.(['incorrect_stim' num2str(j)]),'linewidth',2},0.1);
        end
    end
    
    xlim([0 opt.sta_pre_frames+opt.sta_post_frames+1])
    axis square
    
    text(1,1,['Cell ' num2str(i) ' (ROI ' num2str(cell_struct(i).cnm_idx) ')'],'units','normalized','color','black','Horizontalalignment','right','VerticalAlignment','top')

    if(~isempty(find(opsin_positive_idx==i)))
        text(1,1,['Cell ' num2str(i)],'units','normalized','color','r','Horizontalalignment','right','VerticalAlignment','top')
    end
    box off
    
    % mark tuned cells
    if( cell_struct(i).is_tuned)
        box on
        this_color = trial_color.(['correct_stim' num2str(cell_struct(i).pref)]);
        set(gca,'XColor',this_color,'YColor',this_color,'linewidth',3)
    end
    
    % mark time
    plot([1,1].*opt.sta_gocue_frame, ylim,'color','black','linestyle',':','linewidth',2)
    if ~ opt.flag_use_peak
        plot([avg_frame_range(1),avg_frame_range(end)],[0,0],'color',[.5 .5 .5],'linewidth',2)
    else
        plot([peak_frame_range(1),peak_frame_range(end)],[0,0],'color',[.5 .5 .5],'linewidth',2)
    end
    
    % show auc
    text(0.05,.8,['trial auc ' num2str(cell_struct(i).correct_stimAUC,'%10.2f')],'units','normalized', 'horizontalalignment','left', 'color','black')
    text(0.05,.7,['zscore auc '  num2str(cell_struct(i).correct_stimAUC_zscore,'%10.2f') ],'units','normalized', 'horizontalalignment','left', 'color','black')
    
    % modify x ticks
    xaxisvalues = [0:30:opt.sta_pre_frames+opt.sta_post_frames];
    xticks(xaxisvalues)
    xticklabels(arrayfun(@(x){num2str(x)},(xaxisvalues-opt.sta_trialon_frame)./opt.frame_rate))
    
end
export_fig([fig_save_path filesep plot_cell_type 'ShadedSTATrace_' strrep(caiman_file,'.mat','.png')])


%%  == DIMENSIONALITY REDUCTION ==
%% params for factor analysis (also used for stim_opt so run this)
fa_opt.bin_size = 1;
fa_opt.gocue_bin = floor(opt.sta_gocue_frame/fa_opt.bin_size);
fa_opt.stim_bin = ceil(opt.sta_stim_frame/fa_opt.bin_size);
fa_opt.frame_rate = 30/fa_opt.bin_size;
fa_opt.Fs = opt.frame_rate;
fa_opt.trial_length = opt.trial_length/fa_opt.bin_size;
fa_opt.trial_color = trial_color;
fa_opt.idx_fields = {opt.trigger_idx_fd};

% fa_opt.avg_frames = fa_opt.stim_bin+30:1:fa_opt.gocue_bin;
fa_opt.avg_frames  = [];
easy_trial_idx = 1:10;

%% dpca
% shuffle takes time - try parfor?
if strcmp(opt.method,'dpca')
dpca_opt = opt;
dpca_opt.trial_color = trial_color;
fd_names ={'stim_1_correct','stim_1_incorrect','stim_2_correct','stim_2_incorrect'};

[dpca_struct,traces_idx_struct] = get_dpca_traj_brief(raw_cell_struct,cell_idx_struct.(opt.trigger_idx_fd),fd_names,opt);

plot_pop_vectors(dpca_struct.traj_struct,fd_names,5,dpca_opt,...
    'plot_ylabel','PC level','plot_num_cols',2);
end
%% run factor analysis
if strcmp(opt.method,'fa')
% if FLAG_GOT_CHOICE_AUC
    fa_opt.fd_names ={'stim_1_correct','stim_1_incorrect','stim_2_correct','stim_2_incorrect'};
% else
%     fa_opt.fd_names ={'stim_1_correct','stim_2_correct'};
%     disp('using correct trials only')
% end
fa_opt.plot_fds = fa_opt.fd_names;
fa_opt.m = 3;
fa_opt.IF_MEDFILT = 0;
disp('Rnning factor analysis...')
tic
% excluding first 10 trials that were used for baselining - didnt make a
% difference..
[~,trial_idx_struct ]= cellfun(@(x)setdiff(trial_indices.(x),easy_trial_idx),fa_opt.fd_names,'UniformOutput',false);

[traces_in,fa_trial_idx,num_trials] = get_input_seq(raw_cell_struct,cell_idx_struct.(fa_opt.idx_fields{1}),...
    fa_opt.fd_names,fa_opt.bin_size,'IF_MEDFILT',fa_opt.IF_MEDFILT,...
    'fd_trial_idx',trial_idx_struct,'avg_frames',fa_opt.avg_frames);%% 
fa_struct = struct();
fa_struct.mean = mean(traces_in);
fa_struct.std = std(traces_in);
[fa_struct.lambda,fa_struct.psi,fa_struct.T,fa_struct.stats,fa_struct.F] = factoran(traces_in,fa_opt.m,'Xtype','data','Maxit',1000);
invsqrtPsi = diag(1 ./  sqrt(fa_struct.psi)); % get transition matrix (multiply to zscored data)
fa_struct.transmat = invsqrtPsi/(fa_struct.lambda'*invsqrtPsi);
if isempty(fa_opt.avg_frames)
    [fa_traj_struct,num_fs] = get_pop_vectors(fa_struct.F,fa_opt.trial_length,fa_trial_idx); % ran FA on full traces
else
    full_traces_in  = get_input_seq(raw_cell_struct,cell_idx_struct.(fa_opt.idx_fields{1}),...
        fa_opt.fd_names,fa_opt.bin_size,'IF_MEDFILT',fa_opt.IF_MEDFILT,...
        'fd_trial_idx',trial_idx_struct);%%
    F = (full_traces_in-fa_struct.mean)./fa_struct.std;
    [fa_traj_struct,num_fs] = get_pop_vectors(F,fa_opt.trial_length,fa_trial_idx);
end
plot_pop_vectors(fa_traj_struct,fa_opt.plot_fds,fa_opt.m,fa_opt,...
    'plot_ylabel','Factor level','plot_num_cols',2);
toc
disp('...Done')
end
    
%% coding direction
% using raw onine traces to mimic online condition
if strcmp(opt.method,'cd')
choice_fds = {'stim_1_correct','stim_2_correct'};
cd_opt = opt;
cd_opt.trial_color = trial_color;
coding_direct = cell2mat(arrayfun(@(x)nanmean(mean(raw_cell_struct(x).(choice_fds{1})(opt.sta_avg_frames,:),2)),cell_idx_struct.(opt.trigger_idx_fd),'UniformOutput', false))...
    - cell2mat(arrayfun(@(x)nanmean(mean(raw_cell_struct(x).(choice_fds{2})(opt.sta_avg_frames,:),2)),cell_idx_struct.(opt.trigger_idx_fd),'UniformOutput', false));
% project to coding direction
cd_traj_struct = get_cd_projection(raw_cell_struct(cell_idx_struct.(opt.trigger_idx_fd)),coding_direct,choice_fds);
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
%% get stim decoder 
stim_opt = fa_opt;
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
export_fig([fig_save_path filesep 'StimDecoderPerform_' strrep(caiman_file,'.mat','.png')])

%% get choice decoder (directly from traj, using threshold trials) 
choice_opt = fa_opt;
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
export_fig([fig_save_path filesep 'ChoiceDecoderPerform_' strrep(caiman_file,'.mat','.png')])

%% get projections on stim decoder
stim_thresh = -stim_struct.B(:,1);

switch opt.method
    case 'fa'
        stim_weights = fa_struct.transmat* stim_struct.B(:,2:end)';
        [stim_norm_weights,stim_norm_thresh] = get_norm_weights(stim_weights,stim_thresh,fa_struct.mean,fa_struct.std); % called in CompareDecoderWithAnimal
        %     stim_norm_weights = stim_weights;
        %     stim_norm_thresh = stim_thresh;
    case 'dpca'
        stim_weights = dpca_struct.W* stim_struct.B(:,2:end)';
        [stim_norm_weights,stim_norm_thresh] = get_norm_weights(stim_weights,stim_thresh,dpca_struct.mean,ones(size(dpca_struct.mean))); % called in CompareDecoderWithAnimal
        
    case 'cd'
        stim_weights =coding_direct'* stim_struct.B(2:end)';
        stim_norm_weights = stim_weights;
        stim_norm_thresh = stim_thresh;
end

test_fd_names = {'stim_1_correct','stim_1_incorrect',...
                 'stim_2_correct','stim_2_incorrect',...
};

stim_proj_struct = struct();
[stim_proj_struct] = get_projections(raw_cell_struct(cell_idx_struct.(opt.trigger_idx_fd)),stim_norm_weights,test_fd_names,'proj_struct',stim_proj_struct,'bias',-stim_norm_thresh,'IS_CELL_STRUCT',1);

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
export_fig([fig_save_path filesep 'StimDecoderProj_' strrep(caiman_file,'.mat','.png')])

%% get projections on choice decoder
choice_thresh = -choice_struct.B(:,1);
switch opt.method
    case 'fa'
        choice_weights = fa_struct.transmat* choice_struct.B(:,2:end)';
        [choice_norm_weights,choice_norm_thresh] = get_norm_weights(choice_weights,choice_thresh,fa_struct.mean,fa_struct.std); % called in CompareDecoderWithAnimal
    case 'dpca'
        choice_weights = dpca_struct.W* choice_struct.B(:,2:end)';
        [choice_norm_weights,choice_norm_thresh] = get_norm_weights(choice_weights,choice_thresh,dpca_struct.mean,ones(size(dpca_struct.mean))); % called in CompareDecoderWithAnimal
        
    case 'cd'
        choice_norm_weights = coding_direct'* choice_struct.B(:,2:end)';
        choice_norm_thresh = choice_thresh;
end

test_fd_names = {'stim_1_correct','stim_1_incorrect',...
                 'stim_2_correct','stim_2_incorrect',...
};
choice_proj_struct = struct();
[choice_proj_struct] = get_projections(raw_cell_struct(cell_idx_struct.(opt.trigger_idx_fd)),choice_norm_weights,...
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
export_fig([fig_save_path filesep 'ChoiceDecoderProj_' strrep(caiman_file,'.mat','.png')])


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

%% SELECT CONDITION STIM TYPES  CHANGE THIS
disp('ENTER HERE!')
opt.decoder_to_use = 'stim'; % or choice

condition_pybehav_stimtypes = [1 2]; % stim types in pybehav that will be monitored for photostim;
condition_pybehav_stimvars = [1 1]; % stim types in pybehav that will be monitored for photostim;
condition_pyrtaoi_stimtypes = [1 2]; % stim type 1 will be stimuated when BELOW, stim type 2 will be stimulated when ABOVE
condition_pyrtaoi_stimvars  = [1 1]; % type and var need to match - vars are not used for texture exp, set all to 1

control_pybehav_stimtypes = [1 2]; % stim types in pybehav that will be monitored for photostim but not giving power
control_pybehav_stimvars = [2 2]; % 
control_pyrtaoi_stimtypes = [1 2]; % NO-Power stim: stim type 1 will be stimuated when BELOW, stim type 2 will be stimulated when ABOVE
control_pyrtaoi_stimvars  = [2 2]; % type and var need to match - vars are not used for texture exp, set all to 1



catch_pybehav_stimtypes = [3,4,3,4,3,4]; % 4 with reward, 3 no-reward
catch_pybehav_stimvars = [1,1,2,2,3,3]; % photostim targets. 1:tex1, 2:tex2, 3:none
catch_pyrtaoi_stimtypes = [3,3,4,4,5,5]; % 3 stim tex1, 4 stim tex2, 5 no stim
catch_pyrtaoi_stimvars = [1,1,1,1,1,1];

condition_type = {[101],[201],[301],[401],[102],[202]}; % match target ensembles 100*stim_type + stim_va
opt.target_idx_fd = {'all_tex1','all_tex2','all_tex1','all_tex2','all_tex1','all_tex2'}; % target groups corresponding to condition_type

% opt.target_idx_fd = {'tex1','tex2','tex1','tex2','tex1','tex2'}; % target groups corresponding to condition_type
% above corresponds to: stimulate texture 2 group if population activity is
% going above threshold (towards texutre 1 direction) in catch trials where
% port2 will be rewarded.
% control: if not doing closed loop but just randomly giving reward on
% catch trials, then performance in the normal texutre trials will drop
% (more)

keyboard
switch opt.decoder_to_use
    case 'stim'
        decod_struct = stim_struct; % choose which decoder to use
        pop_weights = stim_norm_weights; % select which weights to use
        pop_thresh = stim_norm_thresh;  % select which thresh to use
        opt.pop_opt = stim_opt;

    case 'choice'
        decod_struct = choice_struct; % choose which decoder to use
        pop_weights = choice_norm_weights; % select which weights to use
        pop_thresh = choice_norm_thresh;  % select which thresh to use
        opt.pop_opt = choice_opt;

end

% generate parameters for pyRTAOI population analysis
pop_params = struct();
pop_params.weights = pop_weights;
pop_params.thresh = pop_thresh;
pop_params.condition_type = condition_type; 
% frames_of_interest = [pop_params.frames_enable_trigger(1):pop_params.frames_enable_trigger(end)]+opt.sta_baseline_frames;
% [proj_sd,fds_of_interest] = get_proj_noise(choice_proj_struct,condition_stim_type,condition_stim_var,frames_of_interest);
% plot_pop_vectors(choice_proj_struct,fds_of_interest,1,opt,...
%        'noise_thresh',proj_sd,'ylimit',[-100 100],'plot_ylabel','proj to choice-stim','plot_num_cols',2,'IF_PLOT_RAW_ONLY',1)
pop_params.thresh_sd = 0;% NOT USING SD
pop_params.frames_enable_trigger  = [100 120]; % monitor frames, frame 120 is go-cue

%% project catch trials on decoder
test_fd_names = {'stim_5_miss','stim_5_lick'};
[catch_proj_struct] = get_projections(raw_cell_struct(cell_idx_struct.(opt.trigger_idx_fd)),pop_weights,...
    test_fd_names,'bias',-pop_thresh,'IS_CELL_STRUCT',1);
% compare projection on choice axis trials
plot_pop_vectors(catch_proj_struct,test_fd_names,1,opt,...
    'plot_ylabel','Projection','plot_num_cols',2,'IF_PLOT_RAW_ONLY',1)
suptitle('Catch trials decoder projections')

%% GET PHOTOEXCITABLE TARGETS - optional
% load result file from OnlinePhotoexcitability.m
% run the script in another matlab!
IF_PHOTOTEST_LOADED = false;
try
    [photo_file,photo_path] = uigetfile('*.mat','Select ProcPhotoExci data',caiman_path);
    disp(['Loaded file :',fullfile(photo_path,photo_file)])
    load(fullfile(photo_path,photo_file));
    IF_PHOTOTEST_LOADED = true;
catch
    warning('photo file not loaded')
end

if IF_PHOTOTEST_LOADED
    % merge cell struct into one
    photo_fds = fields(photo_output_struct.cell_struct);
    for c = 1:num_cells
        for f = 1:numel(photo_fds)
            cell_struct(c).(photo_fds{f}) = photo_output_struct.cell_struct(c).(photo_fds{f});
        end  
    end
    % get photoexcitable target indices
    photo_idx.all = find(extractfield(cell_struct,'is_photo')>0);
    cell_idx_struct.photo_stim1 = intersect(photo_idx.all, cell_idx_struct.(opt.target_idx_fd{1}));
    cell_idx_struct.photo_stim2 = intersect(photo_idx.all, cell_idx_struct.(opt.target_idx_fd{2}));
    opt.target_idx_fd = {'photo_stim1','photo_stim2','photo_stim1','photo_stim2','photo_stim1','photo_stim2'}; % overide target indices with the excitable ones
end
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




%% =================== MAKE OUTPUT FILE FOR PYRTAOI =======================
[output2,save_time2] = generate_cell_idx_file(cell_struct,cell_idx_struct,pop_params,opt);

%% ==================== MAKE PYBEHAV EXTERNAL FILE ========================
% make file for pybehavior
num_condition_per_loop = 5; 
num_control_per_loop = 3;
loop_pybehav_stim_types = [repmat(condition_pybehav_stimtypes,[1,num_condition_per_loop]),...
    repmat(control_pybehav_stimtypes,[1,num_control_per_loop]),catch_pybehav_stimtypes]; % condition trial types plus two extreme stim types
loop_pybehav_stim_vars  = [repmat(condition_pybehav_stimvars,[1,num_condition_per_loop]),...
    repmat(control_pybehav_stimvars,[1,num_control_per_loop]),catch_pybehav_stimvars];
loop_pyrtaoi_stim_types = [repmat(condition_pyrtaoi_stimtypes ,[1,num_condition_per_loop]),...
    repmat(control_pyrtaoi_stimtypes,[1,num_control_per_loop]),catch_pyrtaoi_stimtypes]; % condition trial types plus two extreme stim types
loop_pyrtaoi_stim_vars  = [repmat(condition_pyrtaoi_stimvars ,[1,num_condition_per_loop]),...
    repmat(control_pyrtaoi_stimvars,[1,num_control_per_loop]),catch_pyrtaoi_stimvars];

num_loops = 20; % will loop through loop_stim_types and loop_stim_vars, randomise within each loop
file_save_name = [opt.exp_name '_PyBehavior_' save_time2];
[pybehav_seq,pyrtaoi_seq] = generate_texture_trial_seq(loop_pybehav_stim_types,loop_pybehav_stim_vars,...
    loop_pyrtaoi_stim_types,loop_pyrtaoi_stim_vars,num_loops,opt.output_path,file_save_name);
% make file for pyrtaoi (use for 'load pybahave file')
pyrtaoi_file_save_name = [opt.exp_name '_RTAOiPyBehavior_' save_time2];
dlmwrite([opt.output_path filesep pyrtaoi_file_save_name '.txt'],pyrtaoi_seq)
disp(['saved as:' opt.output_path filesep pyrtaoi_file_save_name '.txt'] )

%% save cell structure 
save_brief_fds = {'stim_1_correct','stim_2_correct','stim_1_incorrect','stim_2_incorrect',...
    'correct_stimAUC','correct_stimAUC_zscore','is_tuned','photo_auc_zscore','is_photo',...
    'cnm_idx','jsf','cnn_prediction'};
save_cell_struct_brief = struct();
for c = 1:num_cells
    for f = 1:numel(save_brief_fds)
        try
        save_cell_struct_brief(c).(save_brief_fds{f}) =  cell_struct(c).(save_brief_fds{f});
        catch
            disp([save_brief_fds{f} 'not saved'])
        end
    end
end

tex_output_struct = struct();
tex_output_struct.pop_params = pop_params;
tex_output_struct.cell_struct = save_cell_struct_brief;
tex_output_struct.cell_idx_struct = cell_idx_struct;
tex_output_struct.opt = opt;
tex_output_struct.input_caiman_file = fullfile(caiman_path,caiman_file);
tex_output_struct.pybehav_seq = pybehav_seq;
tex_output_struct.pyrtaoi_seq = pyrtaoi_seq;


output_save_name = [save_path filesep  'ProcTex_' save_time2 '_' caiman_file ];
save(output_save_name,'tex_output_struct')
disp(['Output struct saved as:' output_save_name ])

%% for debug
trigger_idx = cell_idx_struct.(opt.trigger_idx_fd);
online_sd = caiman_data.sd_level(trigger_idx);
online_bs = caiman_data.bs_level(trigger_idx);
online_w = caiman_data.ROIw(trigger_idx);

%% load struct

%% end for debug

%% =========================== CHECK PLOTS ======================================

%% plot the entire trajectory - check drift of baseline between initialisation movie and condition movie

% temp_filtC = cell2mat({cell_struct(:).('filtC')}');
% online_filtC = nan(num_cells,tot_frames);
% for c = 1:num_cells
% temp_trace = nan(1,tot_frames);
% temp_trace(caiman_frames) = temp_filtC(c,1:num_frames);
% online_filtC(c,:) = fillmissing(temp_trace,'linear',1);
% end
% test_traj = pop_params.weights'*online_filtC(cell_idx_struct.(opt.trigger_idx_fd),:)-pop_params.thresh;
% 
% figure;hold on;
% plot(test_traj);
%%
% %% Plot STAs for trigger cells (cells to monitor)
% figure('name','trigger sta traces')
% num_plot_cols = 4;
% num_plot_rois = length(cell_idx_struct.(opt.trigger_idx_fd));
% num_plot_rows = ceil(num_plot_rois/num_plot_cols);
% plot_count = 1;
% for ii = 1:num_plot_rois
%     i = cell_idx_struct.(opt.trigger_idx_fd)(ii);
%     subtightplot(num_plot_rows,num_plot_cols,plot_count)
%     % plot traces
%     hold on
%     for t = 1:size(cell_struct(i).sta_traces,1)
%         plot(cell_struct(i).sta_traces(t,:),'color',opt.trial_color(t,:) ,'linewidth',1)
%     end
%     plot(cell_struct(i).sta_trace,'color',[.5 .5 .5],'linewidth',1.5)
%     
%     xlim([0 opt.sta_pre_frames+opt.sta_post_frames+1])
%     set(gca,'xtick',[],'xcolor',[1 1 1])
%     axis square
%     
%     plot([opt.sta_pre_frames opt.sta_pre_frames],ylim,'color',[.5 .5 .5]) % start of withold window
%     plot([opt.sta_pre_frames opt.sta_pre_frames] + length(opt.withold_frames_adj),ylim,'color',[0 0 0]) % go-cue
% 
%     plot_count = plot_count+1;
%     
%     text(0.05,1,['ROI ' num2str(i)],'units','normalized', 'horizontalalignment','left', 'color','black')
%     text(0.05,.9,['sensory auc ' num2str(cell_struct(i).correct_stimAUC,'%10.2f')],'units','normalized', 'horizontalalignment','left', 'color','black')
%     text(0.05,.8,['zscore auc '  num2str(cell_struct(i).correct_stimAUC_zscore,'%10.2f') ],'units','normalized', 'horizontalalignment','left', 'color','black')
%     
% end
% suptitle([opt.trigger_idx_fd ' rois, stim-triggered response'])


% %% Show sensory cells on maps
% figure('name','pref. texture on fov');
% subplot(1,2,1)
% imagesc(com_fov)
% colormap(gray)
% colorbar('location','southoutside');
% axis square
% title('Detected ROIs')
% 
% ax1 = subplot(1,2,2)
% value_field = 'pref_tex';
% plot_value_in_rois( cell_struct, value_field,[256 256],ax1,'colorlut',[[1,1,1];opt.type_color],'IF_NORM_PIX',0,'IF_CONTOUR',1,'IF_SHOW_OPSIN',1,'zlimit',[0 4]);
% set(gca,'Ydir','reverse')
% title('Sensory cells (colored by pref. texture)')
% 
% %% Plot auc distributions
% figure('name','correct stim auc')
% subplot(1,2,1)
% hold on
% histogram(extractfield(cell_struct,'correct_stimAUC'),'facecolor',[.7 .7 .7],'edgecolor','none','BinWidth',.05)
% histogram(extractfield(cell_struct(extractfield(cell_struct,'is_tuned')==1),'correct_stimAUC'),'facecolor','none','edgecolor',[0,0,1],'BinWidth',.05)
% xlabel('Tex response auc')
% axis square
% 
% subplot(1,2,2)
% hold on
% histogram(extractfield(cell_struct,'correct_stimAUC_zscore'),'facecolor',[.7 .7 .7],'edgecolor','none','BinWidth',.5)
% histogram(extractfield(cell_struct(extractfield(cell_struct,'is_tuned')==1),'correct_stimAUC_zscore'),'facecolor','none','edgecolor',[0,0,1],'BinWidth',.5)
% xlabel('Tex response auc (zscore)')
% axis square
% 
% 

