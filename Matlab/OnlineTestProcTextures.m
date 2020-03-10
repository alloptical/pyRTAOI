% analyse conditioning session in delection task

%% add path - change this for rig
clear all
close all
clc

% ZZ PC
% addpath(genpath('Y:\zzhang\Python\pyRTAOI\Matlab'));
% cd('Y:\zzhang\Python\pyRTAOI\Matlab');

% BRUKER1
% matlab_set_paths_zz

%% stim parameter - CHANGE THIS
crop_num_trials = 156; % specify number of trials recorded if aborted half way
IF_GO_NOGO = true;
IF_USE_PYRTAOI_STIMTYPE = true; % for condition session this will be different from pybehav

opt.go_stim_types = [2,3]; % stim types with reward channel == 0 in pybehav
pybehav_condition_types = [1,2,3,4]; % for checking files only. pyrtaoi will use [1,2] for closed-loop stim_type; 3,4 for photstim tex1 and tex2 ensembles in catch trials, and 5 for catch trial without photostim

easy_trial_idx = 1:10;


%% params
opt.N = 1.5; % threshold for significant auc
opt.sta_pre_frames = 150; 
opt.sta_post_frames = 120;
opt.sta_baseline_frames = 30; % relative to beginning of sta traces
opt.trial_length = 1+opt.sta_pre_frames+opt.sta_post_frames;

opt.sta_stim_frame = 90; % roughly when first contact occurs
opt.gocue_frame = 120; % relative to trial start
opt.stimon_frame = 90; % relative to trial start
opt.end_rw_frame = 180; % end of response window

% frame indices relative to sta trace
opt.sta_gocue_frame = opt.sta_pre_frames;
opt.sta_trialon_frame = opt.sta_pre_frames-opt.gocue_frame;
opt.sta_avg_frames = [-30:1:0]+opt.sta_gocue_frame; % 1 sec before go-cue
opt.sta_peak_search_range =  [-60:1:0]+opt.sta_gocue_frame;

opt.withold_frames_adj = [60 120]; 

opt.offcell_avg_frames = 30;
opt.sta_amp_thresh = 1;
opt.frame_rate = 30;

opt.flag_use_peak = true; % if false then will use average to get roc

opt.correct_trial_only = false; 

% select which online trace to use
opt.IF_USE_FILTC = true;

% select which dim reduction to use
% init color
[trial_color] = online_tex_init_color();

%% load CAIMAN data
[caiman_file,caiman_path] = uigetfile('*.mat','Select caiman data');
caiman_data = load(fullfile(caiman_path,caiman_file));
disp(['Loaded file :',fullfile(caiman_path,caiman_file)])
opt.ds_factor = caiman_data.ds_factor;
try
opt.photo_enable_frame = double(caiman_data.offsetFrames);
opt.sta_stimon_frame = opt.photo_enable_frame+ opt.sta_baseline_frames;

catch
    opt.photo_enable_frame = [];
    opt.sta_stimon_frame = [];
    warning('no photostim for the session loaded')
end

%% load Pybehavior data
try
    [pb_file,pb_path] = uigetfile('*.mat','Select pybehavior data',caiman_path);
    behavior_data =  load(fullfile(pb_path,pb_file)) ;
    disp(['Loaded file :',fullfile(pb_path,pb_file)])
    
    FLAG_PYBEHAV_LOADED = true;
catch
    FLAG_PYBEHAV_LOADED = false;
end


%% load config file: baseline output struct (file named as 'OutputParams')
[baseline_file,baseline_path] = uigetfile('*.mat','Select baseline OutputParams',caiman_path);
baseline_output = load(fullfile(baseline_path,baseline_file));
disp(['Loaded file :',fullfile(baseline_path,baseline_file)])
decod_struct = baseline_output.output;
norm_weights = decod_struct.trigger_weights;
norm_thresh = decod_struct.trigger_thresh;
trigger_idx = decod_struct.trigger_idx;
target_idx = decod_struct.target_idx;
target_ensembles = decod_struct.target_ensembles;
thresh_sd = decod_struct.thresh_sd;
%% config save path
save_path = [caiman_path filesep 'analysis_files'];
if ~exist(save_path, 'dir')
    mkdir(save_path)
end

fig_save_path = [caiman_path filesep 'figures'];
if ~exist(fig_save_path, 'dir')
    mkdir(fig_save_path)
end
% setup save path and name
opt.output_path = caiman_path;
opt.exp_name = strrep(caiman_file,'.mat','proc');
disp(['analysis files savepath:' save_path])
%% organise data (generate plots for sanity check)
tot_num_trials = min([crop_num_trials,length(caiman_data.trialOrder),numel(behavior_data.results)]);

if(~isempty(caiman_data.stim_frames_caiman))
    sens_stim_frames = caiman_data.sensory_stim_frames+caiman_data.t_init;
else
    sens_stim_frames = [];
end

% discard sens stim frames beyond tottime (deal with senarios when movie is aborted early)
% note that here sen_stim_frames are trial-trigger frames to pybehav
sens_stim_frames(sens_stim_frames>caiman_data.t_cnm-opt.sta_post_frames) = [];
tot_num_trials = min([tot_num_trials,length(sens_stim_frames)]);
sens_stim_frames = sens_stim_frames(1:tot_num_trials);
num_trials = length(sens_stim_frames);

%% check if trial order matches in behavior and caiman file
if FLAG_PYBEHAV_LOADED
    [trials,odd_trial_idx] = make_trials_struct(behavior_data);   
    % discard trials after num_trials   
    trials = structfun(@(x)x(1:num_trials),trials,'UniformOutput',false);
end
% pyrtaoi trial type
trialOrder = caiman_data.trialOrder(1:num_trials); %
trialTypes = unique(trialOrder);
num_stim_type = length(unique(trialOrder)); % orientations/texture
trials.trialOrder = trialOrder;


% only get correct trials
if opt.correct_trial_only && FLAG_PYBEHAV_LOADED
    sens_stim_frames = sens_stim_frames(trials.correct==1);
    trialOrder = trialOrder(trials.correct==1);
elseif ~ FLAG_PYBEHAV_LOADED
    % make a dummy trial struct (setting all trials to correct)
    trials.correct = ones(size(sens_stim_frames));
    trials.incorrect = zeros(size(sens_stim_frames));
    trials.stim_type = trialOrder;
    trials.stim_var = ones(size(trialOrder));
    trials.miss = zeros(size(sens_stim_frames));
end

% roi indices
opsin_positive = caiman_data.opsin_positive;
accepted_idx = caiman_data.accepted_idx+1;
opsin_positive_idx = accepted_idx(opsin_positive>0);


opt.type_color = [trial_color.('stim1');trial_color.('stim2')];
opt.trial_color = trial_color;
%% get trial indices
%% trial types for go-nogo condition
% pyrtaoi: [1 2 3 4 3 4 5 5];
% pybehav: [1 2 3 3 4 4 3 4];
% texture: [1 2 3 3 3 3 3 3]; % by pybehav
% reward:  [0 2 0 0 2 2 2 0]; % by pybehav, 2 is go, 0 is no-go
% target:  [1 2 1 2 1 2 0 0]; % by pyrtaoi
num_trial_types = 8;

[photostim_trial_idx,num_photo_per_trial] = get_trials_with_photostim( caiman_data.sensory_stim_frames, caiman_data.online_photo_frames );
trials.photostim = zeros(1,num_trials);
trials.photostim(photostim_trial_idx(photostim_trial_idx<=num_trials))=1;

trial_indices = struct(); % % get trialtype-outcome indices
all_pyrtaoi_stimtype =  [1 2 3 4 3 4 5 5]; % will be used as variation below
all_pybehav_stimtype =  [1 2 3 3 4 4 3 4];


for v = 1:num_trial_types
    this_var = all_pyrtaoi_stimtype(v);
    this_stim = all_pybehav_stimtype(v);
    trial_indices.(['stim_' num2str(this_stim) '_var_' num2str(this_var) '_correct'  ]) = find(trials.correct==1&trials.stim_type==this_stim&trials.trialOrder == this_var);
    trial_indices.(['stim_' num2str(this_stim) '_var_' num2str(this_var) '_incorrect' ]) = find(trials.incorrect==1&trials.stim_type==this_stim&trials.trialOrder == this_var);
    trial_indices.(['stim_' num2str(this_stim) '_var_' num2str(this_var) '_miss' ]) = find(trials.miss==1&trials.stim_type==this_stim&trials.trialOrder == this_var);
    trial_indices.(['stim_' num2str(this_stim) '_var_' num2str(this_var) '_photostim'  ]) = find(trials.stim_type==this_stim&trials.trialOrder == this_var & trials.photostim ==1);
    trial_indices.(['stim_' num2str(this_stim) '_var_' num2str(this_var) '_nonphotostim'  ]) = find(trials.stim_type==this_stim&trials.trialOrder == this_var & trials.photostim ==0);
  
    if IF_GO_NOGO && ~isempty(intersect(opt.go_stim_types,this_stim))
        trial_indices.(['stim_' num2str(this_stim) '_var_' num2str(this_var) '_incorrect' ]) = trial_indices.(['stim_' num2str(this_stim) '_var_' num2str(this_var) '_miss' ]);
    end
end
% sort catch trials
% using stim_5 as catch; var as target ensemble
trial_indices.('stim_5_var_1') =  find(trials.trialOrder==3&trials.photostim == 1);
trial_indices.('stim_5_var_2') =  find(trials.trialOrder==4&trials.photostim == 1);
trial_indices.('stim_5_var_0') =  find(trials.trialOrder==5&trials.photostim == 0);


% get trials with photostims
if ~isempty( caiman_data.online_photo_frames)
photo_stim_frames =  caiman_data.online_photo_frames + caiman_data.t_init;
photo_stim_frames(photo_stim_frames>caiman_data.t_cnm)=[];
else
    photo_stim_frames = [];
    disp('no photostim found')
end

disp('sorted trial indices')

%% generate config file for control session
% change closed-loop trials to photostim or no photostim trials
% only matters for pyrtaoi
control_trials = trials;

control_trials.trialOrder(control_trials.photostim==1&control_trials.trialOrder<=2) = 2+control_trials.trialOrder(control_trials.photostim==1&control_trials.trialOrder<=2);
control_trials.trialOrder(control_trials.photostim==0) = 5;

% same stim_type and stim_var pairs but shuffled
shuf_idx = randperm(num_trials-numel(easy_trial_idx));
control_idx = [easy_trial_idx,shuf_idx+numel(easy_trial_idx)];
control_trials = structfun(@(x)x(control_idx),control_trials,'UniformOutput',false);
pyrtaoi_seq = [control_trials.trialOrder; ones(1,num_trials)];
pybehav_seq = [control_trials.stim_type; ones(1,num_trials)];


% make trial sequence file for pyrtaoi and pybehav
save_time = datestr(now,'yyyymmdd_HHMM');
pyrtaoi_file_save_name = [opt.exp_name '_RTAOiPyBehavior_CONTROL_' save_time];
dlmwrite([opt.output_path filesep pyrtaoi_file_save_name '.txt'],pyrtaoi_seq)
disp(['saved as:' opt.output_path filesep pyrtaoi_file_save_name '.txt'] )
pyrtaoi_file_save_name = [opt.exp_name '_PyBehavior_CONTROL_' save_time];
dlmwrite([opt.output_path filesep pyrtaoi_file_save_name '.txt'],pybehav_seq)
disp(['saved as:' opt.output_path filesep pyrtaoi_file_save_name '.txt'] )

%% get online trajectory and weights
online_weights = caiman_data.ROIw;
online_thresh = caiman_data.ROIsumThresh;
online_sd = caiman_data.sd_level;
online_bs = caiman_data.bs_level;
online_w = caiman_data.ROIw;

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
    cnm_struct(i).centroid = cm(i,:);
    %     cnm_struct(i).frame_added = find(cnm_struct(i).noisyC >0,1);
    
    % set skipped frames to nan then interpolate
    temp_trace = nan(1,tot_frames);
    temp_trace(caiman_frames) =  caiman_data.online_C(i+1,1:num_frames); % use cnm_C(i,1:num_frames) for offline analysis!
    temp_trace = fillmissing(temp_trace,'linear');
    cnm_struct(i).onlineC_full = temp_trace;
    
    temp_trace = nan(1,tot_frames);
    temp_trace(caiman_frames) =  caiman_data.cnm_C(i,1:num_frames);
    temp_trace = fillmissing(temp_trace,'linear');
    cnm_struct(i).deconvC_full = temp_trace;
    
    temp_trace = nan(1,tot_frames);
    temp_trace(caiman_frames) =  caiman_data.noisyC(i+1,1:num_frames);
    temp_trace = fillmissing(temp_trace,'linear');
    cnm_struct(i).noisy_full = temp_trace;
    % use stimon (start of deflection) frame as 'stim frame' for sta traces
    cnm_struct(i).stim_frames = sens_stim_frames+opt.gocue_frame;
    
end

temp_trace = nan(1,tot_frames);
temp_trace(caiman_frames) =  caiman_data.noisyC(1,1:num_frames);
temp_trace = fillmissing(temp_trace,'linear');
backgroundC = temp_trace;
disp('made cnm_struct')

%% only get accepted cells
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

figure('name','fov','position',[100 100 1200 800])
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
    
    temp_trace = nan(1,tot_frames);
    temp_trace(caiman_frames) =  caiman_data.filt_C(this_idx,1:num_frames);
    temp_trace = fillmissing(temp_trace,'linear');

    cell_struct(i).filtC = temp_trace;

    if(~isempty(find(opsin_positive_idx==i)))
        cell_struct(i).opsin_positive = 1;
    end
end
%% check alignment
figure; hold on
cell_count = 1;
for e = 1:2
    this_idx = target_ensembles{e};
    for i = 1:length(this_idx)
        cell_count = cell_count+1;
    plot(zscore(cell_struct(this_idx(i)).filtC)+cell_count*5,'black');
    end
end



for i = 1:numel(sens_stim_frames)
    this_color = trial_color.(['stim' num2str(trials.stim_type(i) )]);
    plot([sens_stim_frames(i) sens_stim_frames(i)]+opt.stimon_frame,ylim,'color',this_color) % trial-on + roughly first touch
    plot([sens_stim_frames(i) sens_stim_frames(i)]+opt.gocue_frame,ylim,'color',this_color,'linestyle',':') % go-cue
end

for i = 1:numel(photo_stim_frames)
    plot([photo_stim_frames(i) photo_stim_frames(i)],ylim,'color',trial_color.photostim,'linestyle',':') 
end
xlim([caiman_data.t_init tot_frames])

% for i = 1:numel(cnm_struct(74).stim_frames)
%     plot([cnm_struct(74).stim_frames(i) cnm_struct(74).stim_frames(i)],ylim,'color',[0 0 0]) % sta gocue
% end

%     this_online_trace = cell_struct(74).filtC; % online trace med filtered
%     
%     temp_trace = nan(1,tot_frames);
%     temp_trace(caiman_frames) = this_online_trace(1,1:num_frames);
%     this_online_trace = fillmissing(temp_trace,'linear');
%     plot(this_online_trace,'black')
% 

%% plot full traces
figure; hold on
plot_offset = 5;
stim_cell_count = 1;
non_stim_cell_count = 1;

for i = 1:num_cells
    this_cell_trace = zscore(cnm_struct(cell_struct(i).cnm_idx).deconvC_full);
    plot(this_cell_trace+i*plot_offset,'color','black','linewidth',1.5)
    stim_cell_count = stim_cell_count+1;
end
xlabel('Frames')
ylabel('ROI index')
yticks([ 1:num_cells].*plot_offset)
yticklabels(1:num_cells)
xlim([caiman_data.t_init tot_frames])
ylim([0 num_cells].*plot_offset+5)

% background
plot(backgroundC,'color',[.5 .5 .5],'linestyle',':')

for i = 1:numel(sens_stim_frames)
    this_color = trial_color.(['stim' num2str(trials.stim_type(i) )]);
    plot([sens_stim_frames(i) sens_stim_frames(i)]+opt.stimon_frame,ylim,'color',this_color) % trial-on + roughly first touch
    plot([sens_stim_frames(i) sens_stim_frames(i)]+opt.gocue_frame,ylim,'color',this_color,'linestyle',':') % go-cue
end

for i = 1:numel(photo_stim_frames)
    plot([photo_stim_frames(i) photo_stim_frames(i)],ylim,'color',trial_color.photostim) % stim-on
end


%% get stim triggered average STA
for i = 1:num_cells
    this_cell_trace = cnm_struct(cell_struct(i).cnm_idx).deconvC_full;
    this_online_trace = cell_struct(i).filtC; % online trace med filtered
        
    
    this_num_trials = numel(cnm_struct(cell_struct(i).cnm_idx).stim_frames );
    this_sens_stim_frames =  cnm_struct(cell_struct(i).cnm_idx).stim_frames;
    cell_struct(i).is_sensory = 0;
    cell_struct(i).is_offcell = 0;
    cell_struct(i).pref_orient = [];
    cell_struct(i).sta_amp = 0;
    cell_struct(i).sta_traces = [];
    cell_struct(i).sta_trace = [];
    cell_struct(i).accepted = 0;
    
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
disp('got cell_struct sta')

%% Normlalise traces to baseline 
% get baseline and std from the first sets of easy trials
% - looks more consistent with traninig session to normalise this way
% doesnt make sense when sd is close to zero!! - just normalise to baseline
% in case intensity drifts...
cell_bs = cell2mat(arrayfun(@(x)mean(reshape(cell_struct(x).raw_sta_traces(easy_trial_idx,1:opt.sta_baseline_frames),[1,numel(easy_trial_idx)*opt.sta_baseline_frames]),2),1:num_cells,'UniformOutput', false));
cell_sd = cell2mat(arrayfun(@(x)std(reshape(cell_struct(x).raw_sta_traces(easy_trial_idx,:),[1,numel(easy_trial_idx)*opt.trial_length])),1:num_cells,'UniformOutput', false));
cell_mean = cell2mat(arrayfun(@(x)mean(cell_struct(x).raw_sta_trace(:)),1:num_cells,'UniformOutput', false));

for i = 1:num_cells
    cell_struct(i).raw_sta_traces = (cell_struct(i).raw_sta_traces - cell_bs(i));
    cell_struct(i).raw_sta_trace = (cell_struct(i).raw_sta_trace - cell_bs(i));
end
disp('normalised cell_struct sta_traces')

%% STA of online trajectory
online_traj = caiman_data.online_traj;
traj_struct = struct();
temp_trace = nan(1,tot_frames);
temp_trace(caiman_frames) = online_traj(1,1:num_frames);
online_traj = fillmissing(temp_trace,'linear');

[~,~,~,traj_struct.sta_traces,~,~,traj_struct.sta_trace] =...
    make_sta_from_traces(online_traj,this_sens_stim_frames,opt.sta_pre_frames,opt.sta_post_frames,1:opt.sta_baseline_frames);

% sta of backgound component (for alignment check)
[~,~,~,~,~,bg_sta_traces,bg_sta_trace] =...
    make_sta_from_traces(backgroundC,this_sens_stim_frames ,opt.sta_pre_frames,opt.sta_post_frames,1:opt.sta_baseline_frames);
%% plot online trajectory
figure; hold on
plot(online_traj,'color',[.5 .5 .5])

for i = 1:numel(sens_stim_frames)
    this_color = trial_color.(['stim' num2str(trials.stim_type(i) )]);
    plot([sens_stim_frames(i) sens_stim_frames(i)]+opt.stimon_frame,ylim,'color',this_color) % stim-on
    plot([sens_stim_frames(i) sens_stim_frames(i)]+opt.gocue_frame,ylim,'color',this_color,'linestyle',':') % go-cue
end

for i = 1:numel(photo_stim_frames)
    plot([photo_stim_frames(i) photo_stim_frames(i)],ylim,'color',trial_color.photostim,'linestyle',':') % stim-on
end

plot(xlim,[0 0],'color','black')
% xlim([caiman_data.t_init tot_frames])

%% Sort STAs by different trial types
trial_types = fields(trial_indices);
num_sta_traces = size(cell_struct(1).sta_traces,1);
raw_cell_struct = struct();

for i = 1:numel(trial_types)
    this_fd = trial_types{i};
    this_idx = trial_indices.(this_fd);
    this_idx = this_idx(this_idx<min([tot_num_trials,num_sta_traces]));
    this_idx = setdiff(this_idx,easy_trial_idx);
    if ~isempty(this_idx)
        for c = 1:num_cells
            cell_struct(c).(this_fd) = cell_struct(c).sta_traces( this_idx,:)';
            raw_cell_struct(c).(this_fd) = cell_struct(c).raw_sta_traces( this_idx,:)';
            traj_struct.(this_fd) = traj_struct.sta_traces(this_idx,:);
        end      
    end
end
disp('sorted cell_struct sta_traces')



%% Plot STAs traces for certain cell types
figure('name','condition sta traces','units','normalized','outerposition',[0 0 1 1])
stim_fds = {'stim_1_var_1_correct','stim_2_var_2_correct'...
    'stim_1_var_1_photostim','stim_2_var_2_photostim'}; % extreme stim types

peak_frame_range = opt.sta_peak_search_range;
avg_frame_range = opt.sta_avg_frames;


plot_cell_idx = sort(trigger_idx);
plot_num_cells = numel(plot_cell_idx);

num_plot_cols = 8;
num_plot_rows = ceil(plot_num_cells/num_plot_cols);
for ii = 1:plot_num_cells
    subtightplot(num_plot_rows,num_plot_cols,ii)
    i = plot_cell_idx(ii);
    hold on
    % correct trials
    plot(cell_struct(i).(stim_fds{1}),'color',trial_color.stim_1_correct,'linewidth',1)
    plot(cell_struct(i).(stim_fds{2}),'color',trial_color.stim_2_correct,'linewidth',1)
    
    % photostim trials
    plot(cell_struct(i).(stim_fds{3}),':','color',[0 0 0],'linewidth',1)
    plot(cell_struct(i).(stim_fds{4}),':','color',[.5 .5 .5],'linewidth',1)
   
    xlim([0 opt.sta_pre_frames+opt.sta_post_frames+1])
    axis square
    
    text(1,1,['Cell ' num2str(i) ' (ROI ' num2str(cell_struct(i).cnm_idx) ')'],'units','normalized','color','black','Horizontalalignment','right','VerticalAlignment','top')
    text(1,0.9,['W '  num2str(online_w(i),'%0.001f') ],'units','normalized','color','black','Horizontalalignment','right','VerticalAlignment','top')

    box off
    
    
    % mark time
    plot([1,1].*opt.sta_gocue_frame, ylim,'color','black','linestyle',':')
    if ~ opt.flag_use_peak
        plot([avg_frame_range(1),avg_frame_range(end)],[0,0],'color','black','linewidth',2)
    else
        plot([peak_frame_range(1),peak_frame_range(end)],[0,0],'color','black','linewidth',2)
    end
    
%     % mark trigger cell
%     if( any(trigger_idx==i))
%         box on
%         set(gca,'linewidth',3)
%         text(0.05,.8,['weight'  num2str(norm_weights(trigger_idx==i),'%10.2f') ],'units','normalized', 'horizontalalignment','left', 'color','black')
%         
%     end
    
    
    
    % modify x ticks
    xaxisvalues = [0:30:opt.sta_pre_frames+opt.sta_post_frames];
    xticks(xaxisvalues)
    xticklabels(arrayfun(@(x){num2str(x)},(xaxisvalues-opt.sta_trialon_frame)./opt.frame_rate))
    
end
suptitle(strrep(strrep(caiman_file,'.mat',''),'_',' '))
export_fig([fig_save_path filesep 'STATrace_' strrep(caiman_file,'.mat','.png')])
%% get photostim sta amp
photo_ensembles = decod_struct.target_ensembles;
photo_types = floor(cell2mat(decod_struct.condition_type)./100);
catch_photo_types = [3,4]; % catch texture with photostim
catch_nonphoto_type = 5;   % catch texture without photostim
num_photo_ensembles = numel(photo_ensembles);
for i = 1:num_cells
    for e = 1:num_photo_ensembles
        this_photo_type = pybehav_condition_types(photo_types(e));
        cell_struct(i).(['sta_amp_photo_' num2str(photo_types(e))]) = mean(mean(cell_struct(i).sta_traces(trials.photostim==1&trials.trialOrder==this_photo_type,opt.sta_avg_frames)));
        if ~isempty(intersect(this_photo_type,catch_photo_types))
            cell_struct(i).(['sta_amp_nonphoto_' num2str(photo_types(e))]) = mean(mean(cell_struct(i).sta_traces(trials.photostim==0&trials.trialOrder==catch_nonphoto_type,opt.sta_avg_frames)));
        else    
           cell_struct(i).(['sta_amp_nonphoto_' num2str(photo_types(e))]) = mean(mean(cell_struct(i).sta_traces(trials.photostim==0&trials.trialOrder==this_photo_type,opt.sta_avg_frames)));
        end
        cell_struct(i).(['sta_amp_diffphoto_' num2str(photo_types(e))]) = cell_struct(i).(['sta_amp_photo_' num2str(this_photo_type)]) -cell_struct(i).(['sta_amp_nonphoto_' num2str(this_photo_type)]);
    end
    
end
%% Plot photostim STA amp on FOV
for e = 1:num_photo_ensembles
    figure('name',['photstim response on fov type' num2str(e)],'units','normalized','outerposition',[0 0 1 1]);
    plot_count = 1;

    ax = subplot(1,3,plot_count);
    if e==1
        [~,zlimit] = plot_value_in_rois( cell_struct, ['sta_amp_photo_' num2str(photo_types(e))],[256 256],ax,'IF_NORM_PIX',0,...
            'IF_CONTOUR',1,'IF_SHOW_OPSIN',0,'target_cell_idx',photo_ensembles{e},'show_cell_idx',photo_ensembles{e});
        plot_zlimit = zlimit;
    else
        plot_value_in_rois( cell_struct, ['sta_amp_photo_' num2str(photo_types(e))],[256 256],ax,'IF_NORM_PIX',0,...
            'IF_CONTOUR',1,'IF_SHOW_OPSIN',0,'target_cell_idx',photo_ensembles{e},'zlimit',plot_zlimit,'show_cell_idx',photo_ensembles{e});
    end
    title(['Stim type:' num2str(photo_types(e)) ', photo+'])
    plot_count= plot_count+1;
    
    ax = subplot(1,3,plot_count);
    plot_value_in_rois( cell_struct, ['sta_amp_nonphoto_' num2str(photo_types(e))],[256 256],ax,'IF_NORM_PIX',0,...
        'IF_CONTOUR',1,'IF_SHOW_OPSIN',0,'zlimit',plot_zlimit,'target_cell_idx',photo_ensembles{e},'show_cell_idx',photo_ensembles{e});
    plot_count= plot_count+1;
    title(['Stim type:' num2str(photo_types(e)) ', photo-'])
    
    ax = subplot(1,3,plot_count);
    plot_value_in_rois( cell_struct, ['sta_amp_diffphoto_' num2str(photo_types(e))],[256 256],ax,'IF_NORM_PIX',0,...
        'IF_CONTOUR',1,'IF_SHOW_OPSIN',0,'target_cell_idx',photo_ensembles{e},'show_cell_idx',photo_ensembles{e});
    plot_count= plot_count+1;
    title(['Stim type:' num2str(photo_types(e)) ', diff'])
    export_fig([fig_save_path filesep 'PhotoSTA_FOV_Stim' num2str(e) strrep(caiman_file,'.mat','.png')])
    
end

%% Plot photostim STA traces
peak_frame_range = opt.sta_peak_search_range;
avg_frame_range = opt.sta_avg_frames;

plot_cell_idx = target_idx;
plot_num_cells = numel(plot_cell_idx);

num_plot_cols = 8;
num_plot_rows = ceil(plot_num_cells/num_plot_cols);
for photo_idx = 1:num_photo_ensembles
            this_photo_type = pybehav_condition_types(photo_idx);
figure('name','condition sta traces','units','normalized','outerposition',[0 0 1 1])

for ii = 1:plot_num_cells
    subtightplot(num_plot_rows,num_plot_cols,ii)
    i = plot_cell_idx(ii);
    hold on
    % photostimulated trials in black
    plot(cell_struct(i).sta_traces(trials.photostim==1&trials.stim_type==this_photo_type,:)','color',[0 0 0],'linewidth',1)
    plot(cell_struct(i).sta_traces(trials.photostim==0&trials.stim_type==this_photo_type,:)','--','color',[.5 .5 .5],'linewidth',1)
    
    xlim([0 opt.sta_pre_frames+opt.sta_post_frames+1])
    axis square
    
    text(1,1,['Cell ' num2str(i) ' (ROI ' num2str(cell_struct(i).cnm_idx) ')'],'units','normalized','color','black','Horizontalalignment','right','VerticalAlignment','top')
    

    box off
    
    
    % mark time
    plot([1,1].*opt.sta_gocue_frame, ylim,'color','black','linestyle',':')
    if ~ opt.flag_use_peak
        plot([avg_frame_range(1),avg_frame_range(end)],[0,0],'color','black','linewidth',2)
    else
        plot([peak_frame_range(1),peak_frame_range(end)],[0,0],'color','black','linewidth',2)
    end
    
    % mark trigger cell
    if( any(trigger_idx==i))
        box on
        set(gca,'linewidth',3)
        text(0.05,.8,['weight'  num2str(norm_weights(trigger_idx==i),'%10.2f') ],'units','normalized', 'horizontalalignment','left', 'color','black')
    end
    
    % mark target cell in this ensemble
    if( any(photo_ensembles{photo_idx}==i))
        box on
        set(gca,'XColor','r','YColor','r','linewidth',2)
        
    end
    
    % modify x ticks
    xaxisvalues = [0:30:opt.sta_pre_frames+opt.sta_post_frames];
    xticks(xaxisvalues)
    xticklabels(arrayfun(@(x){num2str(x)},(xaxisvalues-opt.sta_trialon_frame)./opt.frame_rate))
    
end
suptitle([strrep(strrep(caiman_file,'.mat',''),'_',' ') ' StimType: ' num2str(photo_idx)])
export_fig([fig_save_path filesep 'Stim' num2str(photo_idx) '_STATrace_' strrep(caiman_file,'.mat','.png')])
end
%% set field names and colors
test_opt = opt;
test_opt.fd_names = {'stim_1_var_1_correct','stim_1_var_1_incorrect',...
    'stim_2_var_2_correct','stim_2_var_2_incorrect'...
};
for i = 1:numel(test_opt.fd_names )
    this_fd = test_opt.fd_names {i};
    if contains(this_fd,'stim_1') && contains(this_fd,'_correct')
        trial_color.(this_fd) = trial_color.correct_stim1;
    end
    if contains(this_fd,'stim_1')&& contains(this_fd,'_incorrect')
        trial_color.(this_fd) = trial_color.incorrect_stim1;
    end
    
    if contains(this_fd,'stim_2') && contains(this_fd,'_correct')
        trial_color.(this_fd) = trial_color.correct_stim2;
    end
    if contains(this_fd,'stim_2')&& contains(this_fd,'_incorrect')
        trial_color.(this_fd) = trial_color.incorrect_stim2;
    end
    
    
    if contains(this_fd,'_miss')
        trial_color.(this_fd) = [94, 34, 92]./255;
    end
    
end
test_opt.trial_color = trial_color;

%% online recorded trajectory
fds_of_interest = {'stim_1_var_1_correct','stim_1_var_1_incorrect','stim_1_var_1_photostim',......
    'stim_2_var_2_correct', 'stim_2_var_2_incorrect','stim_2_var_2_photostim'};
plot_pop_vectors(traj_struct,fds_of_interest,1,test_opt,...
    'noise_thresh',thresh_sd,'plot_ylabel','Projection',...
    'ylimit',[-300 300],...
    'plot_num_cols',3,'IF_PLOT_RAW_ONLY',1)
suptitle('Closed-loop condition trials, online trajectory')

fds_of_interest = {'stim_5_var_1','stim_5_var_2','stim_5_var_0'};  
plot_pop_vectors(traj_struct,fds_of_interest,1,test_opt,...
    'noise_thresh',thresh_sd,'plot_ylabel','Projection',...
    'ylimit',[-300 300],...
    'plot_num_cols',3,'IF_PLOT_RAW_ONLY',1)
suptitle('Catch condition trials, online trajectory')


%% trajectory computed using raw_sta_traces
fds_of_interest = {'stim_1_var_1_correct','stim_1_var_1_incorrect','stim_1_var_1_photostim',......
    'stim_2_var_2_correct', 'stim_2_var_2_incorrect','stim_2_var_2_photostim',...
    'stim_5_var_1','stim_5_var_2','stim_5_var_0'};
[stim_proj_struct] = get_projections(raw_cell_struct(trigger_idx),norm_weights,fds_of_interest,'bias',-norm_thresh,'IS_CELL_STRUCT',1);

plot_pop_vectors(stim_proj_struct,fds_of_interest,1,opt,...
        'ylimit',[-300,300],'plot_ylabel','Projection','plot_num_cols',3,'IF_PLOT_RAW_ONLY',1)
suptitle('Stim decoder projections')
%% plot the entire trajectory - drift of baseline between initialisation movie and condition movie
% looks okay for consecutive movies, but a big jump between movies taken 1.5 hours apart..
temp_filtC = cell2mat({cell_struct(:).('filtC')}');

test_traj = online_w*online_filtC-online_thresh;

figure;hold on;
plot(test_traj);
plot(online_traj);

%%  project to trainer decoder
proj_struct = struct();
[proj_struct] = get_projections(cell_struct(trigger_idx),norm_weights,test_opt.fd_names,'proj_struct',proj_struct,'bias',-norm_thresh,'IS_CELL_STRUCT',1);
plot_pop_vectors(proj_struct,fds_of_interest,1,test_opt,...
    'noise_thresh',thresh_sd,'ylimit',[-5 5],'plot_ylabel','Projection','plot_num_cols',2,'IF_PLOT_RAW_ONLY',1)
suptitle(['Choice decoder projections:' strrep(strrep(caiman_file,'.mat',''),'_',' ')])

%% plot accuracy - select tiral types to condition according to these plots
test_decod_struct = [];
num_compares = round(numel(test_opt.fd_names)/3);
for i = 1:num_compares
    IF_REVERSE = 0;
    this_correct_fd = test_opt.fd_names{3*(i-1)+1};
    this_incorrect_fd = test_opt.fd_names{3*(i-1)+2};
    this_cmp_fds = {this_correct_fd,this_incorrect_fd};
    this_stim_type = strsplit(this_correct_fd,'_');
    this_stim_type = this_stim_type{2};
    % reverse for stim type2
    if strcmp(this_stim_type,'2')||strcmp(this_stim_type,'3')
        IF_REVERSE = 1;
    end
    try
        [  test_decod_struct{i} ] =  get_binary_decoder_disc_time( proj_struct, struct(),...
            this_cmp_fds,test_opt,'IF_FRAMEWISE',0,'threshold',0,'IF_REVERSE',IF_REVERSE);
        test_decod_struct{i}.correct_fd = this_correct_fd;
        test_decod_struct{i}.incorrect_fd = this_incorrect_fd;
        
        figure('name','decoder performance tests')
        plot_binary_decoder(test_decod_struct{i},test_opt)
        suptitle(cellfun(@(x)strrep(x,'_',' '),this_cmp_fds, 'UniformOutput', false))
        
    end
end





%% ============================     END    ================================

