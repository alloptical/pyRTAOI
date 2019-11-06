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
crop_num_trials = 146;
opt.N = 1.5; % threshold for significant auc
opt.sta_pre_frames = 90; % sta will be taken relative to stim-on
opt.sta_post_frames = 120;
opt.sta_baseline_frames = 30; % relative to beginning of sta traces, trial-trigger frame in sta traces
opt.gocue_frame = 105; % relative to trial start, 3.5s
opt.stimon_frame = 60; % relative to trial start
opt.trial_length = 1+opt.sta_pre_frames+opt.sta_post_frames;

% frame indices relative to sta trace
opt.sta_stimon_frame = opt.sta_pre_frames;
opt.sta_trialon_frame = opt.sta_pre_frames-opt.gocue_frame;
opt.sta_avg_frames = [15:1:45]+opt.stimon_frame+opt.sta_baseline_frames;
opt.sta_peak_search_range =  opt.sta_avg_frames;
opt.sta_gocue_frame = opt.gocue_frame+opt.sta_baseline_frames;

opt.gocue_bin = opt.sta_gocue_frame;
opt.stim_bin = opt.sta_stimon_frame;


opt.sta_amp_thresh = 1;
opt.frame_rate = 30;
opt.Fs = opt.frame_rate;

opt.flag_use_peak = false; % if false then will use average to get roc
opt.correct_trial_only = false; % only use correct trials to get tunning

% select cell identity for readout and stimulation
opt.target_idx_fd = {'stim1','stim2'};
opt.trigger_idx_fd = 'port';

opt.plot_stim_types = [1 1 1 1 3 2 2 2 2]; % these pairs should be matched as in experiment_trial_seq
opt.plot_var_types  = [3 5 4 7 2 2 4 5 6];

[trial_color] = deflect_init_color();
%% load CAIMAN data
[caiman_file,caiman_path] = uigetfile('*.mat','Select caiman data');
caiman_data = load(fullfile(caiman_path,caiman_file)); 
disp(['Loaded file :',fullfile(caiman_path,caiman_file)])
opt.ds_factor = caiman_data.ds_factor;
%% load Pybehavior data
try
[pb_file,pb_path] = uigetfile('*.mat','Select pybehavior data',caiman_path);
behavior_data =  load(fullfile(pb_path,pb_file)) ;
disp(['Loaded file :',fullfile(pb_path,pb_file)])
[trials,odd_trial_idx] = make_trials_struct(behavior_data);
FLAG_PYBEHAV_LOADED = true;
catch
    FLAG_PYBEHAV_LOADED = false;
end
%% load trial trigger indices with cue (in case pybehav did not follow trial trigger)
try
[pb_file,pb_path] = uigetfile('*.mat','Select trial trigger idx',caiman_path);
disp(['Loaded file :',fullfile(pb_path,pb_file)])
trialtriggers_idx=  load(fullfile(pb_path,pb_file));
trial_idx = trialtriggers_idx.trialtrigger_idx;
FLAG_TRIALTRIGGER_IDX_LOADED =  true;
tot_num_trials = min([crop_num_trials,length(caiman_data.trialOrder),length(trial_idx),length(trials.stim_type)]); % in case some trial triggers were sent after the session by manual clicking when bebehav is still recording
trial_idx = trial_idx(1:tot_num_trials);
catch
    trial_idx = [];
    tot_num_trials = min([crop_num_trials,length(caiman_data.trialOrder),length(trials.stim_type)]);
    disp('Trusting trial triggers')
    FLAG_TRIALTRIGGER_IDX_LOADED =  false;
end
%% load baseline output struct (file named as 'OutputParams')
[baseline_file,baseline_path] = uigetfile('*.mat','Select baseline OutputParams',caiman_path);
baseline_output = load(fullfile(baseline_path,baseline_file)); 
disp(['Loaded file :',fullfile(baseline_path,baseline_file)])
decod_struct = baseline_output.output;
norm_weights = decod_struct.trigger_weights;
norm_thresh = decod_struct.trigger_thresh;
trigger_idx = decod_struct.trigger_idx;
fds_of_interest = decod_struct.fds_of_interest;
thresh_sd = decod_struct.thresh_sd;
online_weights = caiman_data.ROIw;
online_thresh = caiman_data.ROIsumThresh;
online_sd = caiman_data.sd_level;
online_bs = caiman_data.bs_level;

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

%% check if trial order matches in behavior and caiman file
if FLAG_TRIALTRIGGER_IDX_LOADED % discard caiman trial triggers that didn't evoke a trial
    caiman_data.stim_frames_caiman = caiman_data.stim_frames_caiman(trial_idx);
    caiman_data.sensory_stim_frames = caiman_data.sensory_stim_frames(trial_idx);
    caiman_data.trialOrder = trials.stim_type(1:tot_num_trials); %
    disp('discarded null trial triggers')
end

if FLAG_PYBEHAV_LOADED
    fds = fields(trials);
    for i = 1:numel(fds)
        fd = fds{i};
        trials.(fd) = trials.(fd)(1:tot_num_trials);
    end
    
    if isequal(trials.stim_type,caiman_data.trialOrder)
        disp('data matched')
    else
        warning('trial order mismatch!')
    end
else
    warning('no pybehav data loaded')
end
%% organise data (generate plots for sanity check)
if(~isempty(caiman_data.stim_frames_caiman))
    sens_stim_frames = caiman_data.sensory_stim_frames+caiman_data.t_init;
else
    sens_stim_frames = [];
end

% discard sens stim frames beyond tottime (deal with senarios when movie is aborted early)
% note that here sen_stim_frames are trial-trigger frames to pybehav
sens_stim_frames(sens_stim_frames>caiman_data.t_cnm-opt.sta_post_frames-opt.sta_pre_frames) = [];
num_trials = length(sens_stim_frames);

% discard trials after num_trials
trials = structfun(@(x)x(1:num_trials),trials,'UniformOutput',false);

% trial type
trialOrder = caiman_data.trialOrder(1:num_trials); % this is stim_types, need stim_var to specify deflection amplitudes
trialTypes = unique(trialOrder);
num_stim_type = length(unique(trialOrder)); % orientations/texture

% only get correct trials
if opt.correct_trial_only && FLAG_PYBEHAV_LOADED
    sens_stim_frames = sens_stim_frames(trials.correct==1);
    trialOrder = trialOrder(trials.correct==1);
elseif ~ FLAG_PYBEHAV_LOADED
    % make a dummy trial struct (setting all trials to correct)
    trials.correct = ones(size(sens_stim_frames));
    trials.incorrect = zeros(size(sens_stim_frames));
    trials.stim_type = trialOrder;
end

% roi indices
opsin_positive = caiman_data.opsin_positive;
accepted_idx = caiman_data.accepted_idx+1;
opsin_positive_idx = accepted_idx(opsin_positive>0);

% color trial by stim type
hsv_lut = colormap(hsv);
hsv_lut = hsv_lut(2:end-3,:);
close
indices = round(linspace(1,size(hsv_lut,1),num_stim_type));
opt.trial_color = zeros(numel(trialOrder),3);
for t = 1:numel(trialOrder)
    opt.trial_color(t,:) = trial_color.(['stim' num2str(trialOrder(t))]);
end

opt.type_color = [trial_color.('stim1');trial_color.('stim2')];
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
    
    % set skipped frames to nan
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
    cnm_struct(i).stim_frames = sens_stim_frames+opt.stimon_frame;

end

temp_trace = nan(1,tot_frames);
temp_trace(caiman_frames) =  caiman_data.noisyC(1,1:num_frames);
temp_trace = fillmissing(temp_trace,'linear');
backgroundC = temp_trace;
disp('made cnm_struct')
%% only get accepted cells
accepted_idx = caiman_data.accepted_idx+1;
num_cells = numel(accepted_idx);

%% only analyse frames from the current recording
glob_trialon_frames = caiman_data.sensory_stim_frames + caiman_data.t_init;
glob_trialon_frames(glob_trialon_frames>caiman_data.t_cnm-opt.sta_post_frames-opt.sta_pre_frames) = [];

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
    if(~isempty(find(opsin_positive_idx==i)))
         cell_struct(i).opsin_positive = 1;
    end
end

%% plot full traces
figure; hold on
plot_offset = 5;
stim_cell_count = 1;
non_stim_cell_count = 1;

for i = 1:num_cells
    this_cell_trace = cnm_struct(cell_struct(i).cnm_idx).onlineC_full;
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

for i = 1:numel(glob_trialon_frames)
    this_color = trial_color.(['stim' num2str(trials.stim_type(i) )]);
    plot([glob_trialon_frames(i) glob_trialon_frames(i)]+opt.stimon_frame,ylim,'color',this_color) % trial-on 
    plot([glob_trialon_frames(i) glob_trialon_frames(i)]+opt.gocue_frame,ylim,'color',this_color,'linestyle',':') % go-cue
end

%% get trial indices
% a trial type is defined by stim_type (which lickport is rewarded) and
% stim_var (deflection amplitude)

% get trials with photostims
photo_stim_frames =  caiman_data.online_photo_frames + caiman_data.t_init;
photo_stim_frames(photo_stim_frames>caiman_data.t_cnm)=[];
[photostim_trial_idx,num_photo_per_trial] = get_trials_with_photostim(sens_stim_frames,photo_stim_frames);
trials.photostim = zeros(1,num_trials);
trials.photostim(photostim_trial_idx)=1;

trial_indices = struct(); % % get trialtype-outcome indices
all_trial_types = unique(trials.stim_type.*100+trials.stim_var);
num_trial_types = length(all_trial_types);
all_trial_var = nan(1,num_trial_types);
all_trial_stim =  nan(1,num_trial_types);
for v = 1:num_trial_types
all_trial_var(v) = rem(all_trial_types(v),100);
all_trial_stim(v) = floor((all_trial_types(v)-all_trial_var(v))/100);
end
[all_psycho_struct,num_all_types] = get_psycho_curve(trials,all_trial_stim,all_trial_var);
% [stim var outcome]
for v = 1:num_trial_types
    this_var = all_trial_var(v);
    this_stim = all_trial_stim(v);
    trial_indices.(['stim_' num2str(this_stim) '_var_' num2str(this_var) '_correct'  ]) = find(trials.correct==1&trials.stim_type==this_stim&trials.stim_var == this_var);
    trial_indices.(['stim_' num2str(this_stim) '_var_' num2str(this_var) '_incorrect' ]) = find(trials.incorrect==1&trials.stim_type==this_stim&trials.stim_var == this_var);
    trial_indices.(['stim_' num2str(this_stim) '_var_' num2str(this_var) '_miss' ]) = find(trials.miss==1&trials.stim_type==this_stim&trials.stim_var == this_var);
    trial_indices.(['stim_' num2str(this_stim) '_var_' num2str(this_var) '_photostim'  ]) = intersect(find(trials.stim_type==this_stim&trials.stim_var == this_var),photostim_trial_idx);
end

% photostim trial indices
photo_trials = struct();
fds = fields(trials);
for f = 1:numel(fds)
    fd = fds{f};
    photo_trials.(fd) = double(trials.(fd)).*trials.photostim;
end
[photo_psycho_struct,num_photo_types] = get_psycho_curve(photo_trials,all_trial_stim,all_trial_var);


%% plot psychometric curve
if FLAG_PYBEHAV_LOADED
    
    % [plot_psycho_struct, num_types ]= get_psycho_curve(trials,all_psycho_struct.stim_types,all_psycho_struct.stim_vars);
    plot_psycho_struct = all_psycho_struct;
    num_types = num_all_types;
    
    figure('name','psychometric curve')
    
    f = subplot(3,1,1); hold on
    yyaxis right
    plot(plot_psycho_struct.pc_misses,'-o','color',[.5 .5 .5])
    ylabel('Fraction of miss trials')
    ylim([0 1])
    xlim([0 num_types+1])
    set(gca,'YColor',[0.5 0.5 0.5]);
    
    yyaxis left
    plot(plot_psycho_struct.pc_leftchoices,'-o','color',trial_color.L1)
    plot(plot_psycho_struct.pc_rightchoices,'-o','color',trial_color.L2)
    plot(photo_psycho_struct.pc_leftchoices,'-*','color',trial_color.L1)
    plot(photo_psycho_struct.pc_rightchoices,'-*','color',trial_color.L2)
    
    xlim([0 num_types+1])
    ylim([0 1])
    ylabel('Fraction of choices')
    set(gca,'YColor',[0 0 0]);
    
    subplot(3,1,3); hold on
    b=bar([plot_psycho_struct.stim_types;plot_psycho_struct.stim_vars]');
    b(1).FaceColor=[0 0 0];b(2).FaceColor=[.5 .5 .5];
    legend('Stim type', 'Var type')
    grid on
    ylim([0 9])
    xlim([0 num_types+1])
    xlabel('Trial types')
    
    
    subplot(3,1,2); hold on
    plot(plot_psycho_struct.leftvolts,'-o','color',trial_color.L1)
    plot(plot_psycho_struct.rightvolts,'-o','color',trial_color.L2)
    plot(mean([plot_psycho_struct.rightvolts;plot_psycho_struct.leftvolts],1),'-o','color',[.5 .5 .5])
    xlim([0 num_types+1])
    
    legend('Left','Right')
    ylabel('Deflection strength (V)')
end


%% get stim triggered average 
for i = 1:num_cells
    this_cell_trace = cnm_struct(cell_struct(i).cnm_idx).onlineC_full;
    this_num_trials = numel(cnm_struct(cell_struct(i).cnm_idx).stim_frames );
    this_sens_stim_frames =  cnm_struct(cell_struct(i).cnm_idx).stim_frames;
    cell_struct(i).num_trials = this_num_trials*cell_struct(i).opsin_positive;
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

            cell_struct(i).sta_amp = mean(cell_struct(i).sta_trace(opt.sta_avg_frames));
        
        if  cell_struct(i).sta_amp > opt.sta_amp_thresh
            cell_struct(i).is_sensory = 1;
        end        
    end
    

end
disp('got cell_struct sta')

%% Normlalise traces to baseline and std
% get baseline and std from the first sets of easy trials
% - looks more consistent with traninig session to normalise this way
easy_trial_idx = 1:10;
cell_bs = cell2mat(arrayfun(@(x)mean(reshape(cell_struct(x).sta_traces(easy_trial_idx,1:opt.sta_baseline_frames),[1,numel(easy_trial_idx)*opt.sta_baseline_frames]),2),1:num_cells,'UniformOutput', false));
cell_sd = cell2mat(arrayfun(@(x)std(reshape(cell_struct(x).sta_traces(easy_trial_idx,:),[1,numel(easy_trial_idx)*opt.trial_length])),1:num_cells,'UniformOutput', false));
cell_mean = cell2mat(arrayfun(@(x)mean(cell_struct(x).sta_trace(:)),1:num_cells,'UniformOutput', false));

% using mean and sd calculated online - then got the same value
% cell_sd = online_sd;
% cell_bs = online_bs;
for i = 1:num_cells
    cell_struct(i).sta_traces_zscore = (cell_struct(i).sta_traces - cell_mean(i))./cell_sd(i);
    cell_struct(i).sta_traces = (cell_struct(i).sta_traces - cell_bs(i))./cell_sd(i);
    cell_struct(i).sta_trace = (cell_struct(i).sta_trace - cell_bs(i))./cell_sd(i);
end
disp('normalised cell_struct sta_traces')

%% Sort STAs by different trial types
trial_types = fields(trial_indices);
num_sta_traces = size(cell_struct(1).sta_traces,1);
for i = 1:numel(trial_types)
    this_fd = trial_types{i};
    this_idx = trial_indices.(this_fd);
    this_idx = this_idx(this_idx<min([tot_num_trials,num_sta_traces]));
    if ~isempty(this_idx)
            for c = 1:num_cells
                 cell_struct(c).(this_fd) = cell_struct(c).sta_traces( this_idx,:)';
                 traj_struct.(this_fd) = traj_struct.sta_traces(this_idx,:);
            end
    end
end
disp('sorted cell_struct sta_traces')

%% sta of online trajectory
traj_struct = struct();
online_traj = caiman_data.online_traj;
temp_trace = nan(1,tot_frames);
temp_trace(caiman_frames) = online_traj(1,1:num_frames);
online_traj = fillmissing(temp_trace,'linear');

[~,~,~,~,~,traj_struct.sta_traces,traj_struct.sta_trace] =...
    make_sta_from_traces(online_traj,this_sens_stim_frames,opt.sta_pre_frames,opt.sta_post_frames,1:opt.sta_baseline_frames);

% sta of backgound component (for alignment check)
[~,~,~,~,~,bg_sta_traces,bg_sta_trace] =...
    make_sta_from_traces(backgroundC,this_sens_stim_frames ,opt.sta_pre_frames,opt.sta_post_frames,1:opt.sta_baseline_frames);

%% Plot STAs for all components 
figure('name','condition sta traces','units','normalized','outerposition',[0 0 1 1])
stim_fds = {'stim_1_var_9_correct','stim_2_var_9_correct'}; % extreme stim types
peak_frame_range = opt.sta_peak_search_range;
avg_frame_range = opt.sta_avg_frames;
plot_num_cells = num_cells;
num_plot_cols = 8;
num_plot_rows = ceil(plot_num_cells/num_plot_cols);
for i = 1:plot_num_cells
    subtightplot(num_plot_rows,num_plot_cols,i)
    hold on
    % correct trials

    plot(cell_struct(i).(stim_fds{1}),'color',trial_color.L1,'linewidth',1)
    plot(cell_struct(i).(stim_fds{2}),'color',trial_color.L2,'linewidth',1)

    xlim([0 opt.sta_pre_frames+opt.sta_post_frames+1])
    ylim([-2 10])
    axis square
    
    text(1,1,['Cell ' num2str(i) ' (ROI ' num2str(cell_struct(i).cnm_idx) ')'],'units','normalized','color','black','Horizontalalignment','right','VerticalAlignment','top')

    if(~isempty(find(opsin_positive_idx==i)))
        text(1,1,['Cell ' num2str(i)],'units','normalized','color','r','Horizontalalignment','right','VerticalAlignment','top')
    end
    box off
    
    
    % mark time
    plot([1,1].*opt.sta_stimon_frame, ylim,'color','black','linestyle',':')
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


    
    % modify x ticks
    xaxisvalues = [0:45:opt.sta_pre_frames+opt.sta_post_frames];
    xticks(xaxisvalues)
    xticklabels(arrayfun(@(x){num2str(x)},(xaxisvalues-opt.sta_trialon_frame)./opt.frame_rate))
    
end
suptitle(strrep(strrep(caiman_file,'.mat',''),'_',' '))
export_fig([fig_save_path filesep 'STATrace_' strrep(caiman_file,'.mat','.png')])
%% get photostim sta amp
photo_ensembles = decod_struct.target_ensembles;
photo_vars = rem(cell2mat(decod_struct.condition_type),100);
photo_types = floor(cell2mat(decod_struct.condition_type)./100);

num_photo_ensembles = numel(photo_ensembles);
for i = 1:num_cells
    for e = 1:num_photo_ensembles
        cell_struct(i).(['sta_amp_photo_' num2str(photo_types(e))]) = mean(mean(cell_struct(i).sta_traces(trials.photostim==1&trials.stim_type==photo_types(e),opt.sta_avg_frames)));
        cell_struct(i).(['sta_amp_nonphoto_' num2str(photo_types(e))]) = mean(mean(cell_struct(i).sta_traces(trials.photostim==0&trials.stim_type==photo_types(e),opt.sta_avg_frames)));
        cell_struct(i).(['sta_amp_diffphoto_' num2str(photo_types(e))]) = cell_struct(i).(['sta_amp_photo_' num2str(photo_types(e))]) -cell_struct(i).(['sta_amp_nonphoto_' num2str(photo_types(e))]);
    end

end
%% Plot photostim STA amp on FOV
figure('name','photstim response on fov','units','normalized','outerposition',[0 0 1 1]);
plot_count = 1;
for e = 1:num_photo_ensembles
    ax = subplot(num_photo_ensembles,3,plot_count);
    if e==1
        [~,zlimit] = plot_value_in_rois( cell_struct, ['sta_amp_photo_' num2str(photo_types(e))],[256 256],ax,'IF_NORM_PIX',0,...
        'IF_CONTOUR',0,'IF_SHOW_OPSIN',0,'target_cell_idx',photo_ensembles{e});
        plot_zlimit = zlimit;
    else
        plot_value_in_rois( cell_struct, ['sta_amp_photo_' num2str(photo_types(e))],[256 256],ax,'IF_NORM_PIX',0,...
        'IF_CONTOUR',0,'IF_SHOW_OPSIN',0,'target_cell_idx',photo_ensembles{e},'zlimit',plot_zlimit);
    end
    title(['Stim type:' num2str(photo_types(e)) ', photo+'])
    plot_count= plot_count+1;
    
    ax = subplot(num_photo_ensembles,3,plot_count);
    plot_value_in_rois( cell_struct, ['sta_amp_nonphoto_' num2str(photo_types(e))],[256 256],ax,'IF_NORM_PIX',0,...
        'IF_CONTOUR',0,'IF_SHOW_OPSIN',0,'zlimit',plot_zlimit);
    plot_count= plot_count+1;
    title(['Stim type:' num2str(photo_types(e)) ', photo-'])

    ax = subplot(num_photo_ensembles,3,plot_count);
    plot_value_in_rois( cell_struct, ['sta_amp_diffphoto_' num2str(photo_types(e))],[256 256],ax,'IF_NORM_PIX',0,...
        'IF_CONTOUR',0,'IF_SHOW_OPSIN',0,'target_cell_idx',photo_ensembles{e});
    plot_count= plot_count+1;
    title(['Stim type:' num2str(photo_types(e)) ', diff'])

end
export_fig([fig_save_path filesep 'PhotoSTA_FOV_' strrep(caiman_file,'.mat','.png')])

%% Plot photostim STA traces
figure('name','condition sta traces','units','normalized','outerposition',[0 0 1 1])
peak_frame_range = opt.sta_peak_search_range;
avg_frame_range = opt.sta_avg_frames;
plot_num_cells = num_cells;
num_plot_cols = 8;
num_plot_rows = ceil(plot_num_cells/num_plot_cols);
photo_idx = 1;
for i = 1:plot_num_cells
    subtightplot(num_plot_rows,num_plot_cols,i)
    hold on
    plot(cell_struct(i).sta_traces(trials.photostim==1&trials.stim_type==photo_idx,:)','color',[0 0 0],'linewidth',1)
    plot(cell_struct(i).sta_traces(trials.photostim==0&trials.stim_type==photo_idx,:)','--','color',[.5 .5 .5],'linewidth',1)

    xlim([0 opt.sta_pre_frames+opt.sta_post_frames+1])
    ylim([-2 10])
    axis square
    
    text(1,1,['Cell ' num2str(i) ' (ROI ' num2str(cell_struct(i).cnm_idx) ')'],'units','normalized','color','black','Horizontalalignment','right','VerticalAlignment','top')

    if(~isempty(find(opsin_positive_idx==i)))
        text(1,1,['Cell ' num2str(i)],'units','normalized','color','r','Horizontalalignment','right','VerticalAlignment','top')
    end
    box off
    
    
    % mark time
    plot([1,1].*opt.sta_stimon_frame, ylim,'color','black','linestyle',':')
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
    
    % mark target cell
    if( any(photo_ensembles{photo_idx}==i))
        box on
        set(gca,'XColor','r','YColor','r','linewidth',2)
        
    end

    % modify x ticks
    xaxisvalues = [0:45:opt.sta_pre_frames+opt.sta_post_frames];
    xticks(xaxisvalues)
    xticklabels(arrayfun(@(x){num2str(x)},(xaxisvalues-opt.sta_trialon_frame)./opt.frame_rate))
    
end
suptitle([strrep(strrep(caiman_file,'.mat',''),'_',' ') ' StimType: ' num2str(photo_idx)])
export_fig([fig_save_path filesep 'Stim' num2str(photo_idx) '_STATrace_' strrep(caiman_file,'.mat','.png')])
photo_idx = 2;
for i = 1:plot_num_cells
    subtightplot(num_plot_rows,num_plot_cols,i)
    hold on
    plot(cell_struct(i).sta_traces(trials.photostim==1&trials.stim_type==photo_idx,:)','color',[0 0 0],'linewidth',1)
    plot(cell_struct(i).sta_traces(trials.photostim==0&trials.stim_type==photo_idx,:)','--','color',[.5 .5 .5],'linewidth',1)

    xlim([0 opt.sta_pre_frames+opt.sta_post_frames+1])
    ylim([-2 10])
    axis square
    
    text(1,1,['Cell ' num2str(i) ' (ROI ' num2str(cell_struct(i).cnm_idx) ')'],'units','normalized','color','black','Horizontalalignment','right','VerticalAlignment','top')

    if(~isempty(find(opsin_positive_idx==i)))
        text(1,1,['Cell ' num2str(i)],'units','normalized','color','r','Horizontalalignment','right','VerticalAlignment','top')
    end
    box off
    
    
    % mark time
    plot([1,1].*opt.sta_stimon_frame, ylim,'color','black','linestyle',':')
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
    
    % mark target cell
    if( any(photo_ensembles{photo_idx}==i))
        box on
        set(gca,'XColor','r','YColor','r','linewidth',2)
        
    end

    % modify x ticks
    xaxisvalues = [0:45:opt.sta_pre_frames+opt.sta_post_frames];
    xticks(xaxisvalues)
    xticklabels(arrayfun(@(x){num2str(x)},(xaxisvalues-opt.sta_trialon_frame)./opt.frame_rate))
    
end
suptitle([strrep(strrep(caiman_file,'.mat',''),'_',' ') ' StimType: ' num2str(photo_idx)])
export_fig([fig_save_path filesep 'Stim' num2str(photo_idx) '_STATrace_' strrep(caiman_file,'.mat','.png')])

%% set field names and colors
test_opt = opt;
test_opt.fd_names = {'stim_1_var_3_correct','stim_1_var_3_incorrect','stim_1_var_3_miss',...
                     'stim_1_var_4_correct','stim_1_var_4_incorrect','stim_1_var_4_miss',...
                     'stim_1_var_5_correct','stim_1_var_5_incorrect','stim_1_var_5_miss',...
                     'stim_1_var_7_correct','stim_1_var_7_incorrect','stim_1_var_7_miss',...
                     'stim_3_var_2_correct','stim_3_var_2_incorrect','stim_3_var_2_miss'...
                     'stim_2_var_2_correct','stim_2_var_2_incorrect','stim_2_var_2_miss'...
                     'stim_2_var_4_correct','stim_2_var_4_incorrect','stim_2_var_4_miss'...
                     'stim_2_var_5_correct','stim_2_var_5_incorrect','stim_2_var_5_miss'...
                     'stim_2_var_6_correct','stim_2_var_6_incorrect','stim_2_var_6_miss'};
 for i = 1:numel(test_opt.fd_names )
    this_fd = test_opt.fd_names {i};
    if contains(this_fd,'stim_1') && contains(this_fd,'_correct')
        trial_color.(this_fd) = trial_color.L1;
    end
    if contains(this_fd,'stim_1')&& contains(this_fd,'_incorrect')
        trial_color.(this_fd) = trial_color.ipsi_L1;
    end
    
    if contains(this_fd,'stim_2') && contains(this_fd,'_correct')
        trial_color.(this_fd) = trial_color.L2;
    end
    if contains(this_fd,'stim_2')&& contains(this_fd,'_incorrect')
        trial_color.(this_fd) = trial_color.contra_L2;
    end
    
    if contains(this_fd,'stim_3') && contains(this_fd,'_correct')
        trial_color.(this_fd) = trial_color.correct;
    end
    if contains(this_fd,'stim_3')&& contains(this_fd,'_incorrect')
        trial_color.(this_fd) = trial_color.incorrect;
    end
    
     if contains(this_fd,'stim_1') && contains(this_fd,'_miss')
         trial_color.(this_fd) = [94, 34, 92]./255;
     end
     
     if contains(this_fd,'stim_2') && contains(this_fd,'_miss')
         trial_color.(this_fd) = [34, 64, 94]./255;
     end
     
     if contains(this_fd,'stim_3') && contains(this_fd,'_miss')
         trial_color.(this_fd) = [.3 .3 .3];
     end
end
test_opt.trial_color = trial_color;

%% online recorded trajectory
plot_pop_vectors(traj_struct,fds_of_interest,1,test_opt,...
         'noise_thresh',thresh_sd,'ylimit',[-10 10],'plot_ylabel','Projection','plot_num_cols',2,'IF_PLOT_RAW_ONLY',1)

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

