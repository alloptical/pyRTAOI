% analyse baseline imaging session in delection task
% note: rois idx doesn't match with what's shown in pyRTAOI image window
% run this to quickly find sensory-opsin-positve cells  

% ZZ 2018

% adapted from OnlineProcVisual
% TO DO
% 1. show hit and fa rate for test trials - choice decoder didnt work for
% test trials, try normalise to session mean and std, or to trial baseline
% 2. get noise from trajectory - add that to threshold

% 2. too many manual changes - group them to one
% 3. balance two photostim ensembles 


%% add path - change this for rig
clear all
close all
clc

% ZZ PC
% addpath(genpath('C:\Users\User\Desktop\pyRTAOI-rig\Matlab'));
% cd('C:\Users\User\Desktop\pyRTAOI-rig\Matlab');

% BRUKER1
%  matlab_set_paths_zz

%% parameters - CHANGE THIS
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
opt.sta_avg_frames = [15:1:30]+opt.stimon_frame+opt.sta_baseline_frames;
opt.sta_peak_search_range =  opt.sta_avg_frames;
opt.sta_gocue_frame = opt.gocue_frame+opt.sta_baseline_frames;


opt.sta_amp_thresh = 1;
opt.frame_rate = 30;
opt.Fs = 30;

opt.flag_use_peak = false; % if false then will use average to get roc
opt.correct_trial_only = false; % only use correct trials to get tunning

% select cell identity for readout and stimulation
opt.target_idx_fd = {'stim1','stim2'};
opt.trigger_idx_fd = 'all';

% select population analysis method
opt.IF_USE_HMM = 0;
opt.IF_USE_FA = 1;

% init color
[trial_color] = deflect_init_color();

% select which online trace to use
opt.IF_USE_FILTC = true;

% set minimum snr and cell size - for ROI quality control
opt.min_cell_snr = 0;
opt.min_cell_size = 15; 



%% load CAIMAN data
[caiman_file,caiman_path] = uigetfile('*.mat','Select caiman data');
caiman_data = load(fullfile(caiman_path,caiman_file)); 
disp(['Loaded file :',fullfile(caiman_path,caiman_file)])
opt.ds_factor = caiman_data.ds_factor;
%% load Pybehavior data
try
[pb_file,pb_path] = uigetfile('*.mat','Select pybehavior data',caiman_path);
behavior_data =  load(fullfile(pb_path,pb_file));
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
tot_num_trials = min([length(caiman_data.trialOrder),length(trial_idx),length(trials.stim_type)]); % in case some trial triggers were sent after the session by manual clicking when bebehav is still recording
catch
    trial_idx = [];
    tot_num_trials = min(length(caiman_data.trialOrder),length(trials.stim_type));
    disp('Trusting trial triggers')
    FLAG_TRIALTRIGGER_IDX_LOADED =  false;
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
%% setup save path and name
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
tex_stim_frames = {};
if(~isempty(caiman_data.stim_frames_caiman))
    sens_stim_frames = caiman_data.sensory_stim_frames+caiman_data.t_init;
else
    sens_stim_frames = [];
end

% discard sens stim frames beyond tottime (deal with senarios when movie is aborted early)
sens_stim_frames(sens_stim_frames>caiman_data.t_cnm-opt.sta_post_frames-opt.sta_pre_frames) = [];
num_trials = length(sens_stim_frames);

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

for i = 1:numel(trialTypes)
    tex_stim_frames{i} = sens_stim_frames(trialOrder == trialTypes(i));
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
    cnm_struct(i).onlineC = caiman_data.online_C(i+1,1:num_frames);
    cnm_struct(i).deconvC = caiman_data.cnm_C(i,1:num_frames);
    cnm_struct(i).centroid = cm(i,:);
    
    % set skipped frames to nan
    temp_trace = nan(1,tot_frames);
    temp_trace(caiman_frames) =  caiman_data.cnm_C(i,1:num_frames);
    temp_trace = fillmissing(temp_trace,'linear');
    cnm_struct(i).deconvC_full = temp_trace;
    
    temp_trace = nan(1,tot_frames);
    temp_trace(caiman_frames) =  caiman_data.online_C(i+1,1:num_frames);
    temp_trace = fillmissing(temp_trace,'linear');
    cnm_struct(i).onlineC_full = temp_trace;
    
    % use stimon (start of deflection) frame as 'stim frame' for sta traces
    cnm_struct(i).stim_frames = sens_stim_frames+opt.stimon_frame;

end

temp_trace = nan(1,tot_frames);
temp_trace(caiman_frames) =  caiman_data.online_C(1,1:num_frames);
temp_trace = fillmissing(temp_trace,'linear');
backgroundC = temp_trace;

% only get accepted cells
accepted_idx = caiman_data.accepted_idx+1;
num_cells = numel(accepted_idx);

% only analyse frames from the current recording
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
    
    cell_struct(i).filtC = caiman_data.filt_C(i,1:num_frames);
    
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
    this_cell_trace = cnm_struct(cell_struct(i).cnm_idx).deconvC_full;
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
trial_indices = struct(); % % get trialtype-outcome indices
all_trial_types = unique(trials.stim_type.*100+trials.stim_var);
num_trial_types = length(all_trial_types);
all_trial_var = nan(1,num_trial_types);
all_trial_stim =  nan(1,num_trial_types);
for v = 1:num_trial_types
all_trial_var(v) = rem(all_trial_types(v),100);
all_trial_stim(v) = floor((all_trial_types(v)-all_trial_var(v))/100);
end
all_psycho_struct = get_psycho_curve(trials,all_trial_stim,all_trial_var);
% [stim var outcome]
for v = 1:num_trial_types
    this_var = all_trial_var(v);
    this_stim = all_trial_stim(v);
    trial_indices.(['stim_' num2str(this_stim) '_var_' num2str(this_var) '_correct'  ]) = find(trials.correct==1&trials.stim_type==this_stim&trials.stim_var == this_var);
    trial_indices.(['stim_' num2str(this_stim) '_var_' num2str(this_var) '_incorrect' ]) = find(trials.incorrect==1&trials.stim_type==this_stim&trials.stim_var == this_var);
    trial_indices.(['stim_' num2str(this_stim) '_var_' num2str(this_var) '_miss' ]) = find(trials.miss==1&trials.stim_type==this_stim&trials.stim_var == this_var);

end

%% PSYCHOMETRIC CURVE (if pybehav result given)
if FLAG_PYBEHAV_LOADED
    [plot_psycho_struct, num_types ]= get_psycho_curve(trials,all_psycho_struct.stim_types,all_psycho_struct.stim_vars);
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
    xlim([0 num_types+1])
    ylim([0 1])
    ylabel('Fraction of choices')
    set(gca,'YColor',[0 0 0]);
    
    subplot(3,1,2); hold on
    b=bar([plot_psycho_struct.stim_types;plot_psycho_struct.stim_vars]');
    b(1).FaceColor=[0 0 0];b(2).FaceColor=[.5 .5 .5];
    legend('Stim type', 'Var type')
    grid on
    ylim([0 9])
    xlim([0 num_types+1])
    
    
    subplot(3,1,3); hold on
    plot(plot_psycho_struct.leftvolts,'-o','color',trial_color.L1)
    plot(plot_psycho_struct.rightvolts,'-o','color',trial_color.L2)
    plot(mean([plot_psycho_struct.rightvolts;plot_psycho_struct.leftvolts],1),'--','color',[.5 .5 .5])
    xlim([0 num_types+1])
    
    legend('Left','Right')
    ylabel('Deflection strength (V)')
    xlabel('Trial types')
else
    % plot all trial types and deflection volts
    figure('name','all trial types and voltages');
    subplot(2,1,1)
    hold on
    plot(all_psycho_struct.leftvolts,'-o','color',trial_color.L1)
    plot(all_psycho_struct.rightvolts,'-o','color',trial_color.L2)
    plot(mean([all_psycho_struct.rightvolts;all_psycho_struct.leftvolts],1),'--','color',[.5 .5 .5])
    ylabel('Amplitude (V)')
    legend('Left','Right','avg.')
    xlim([0 num_trial_types+1])
    
    subplot(2,1,2)
    hold on
    b=bar([all_psycho_struct.stim_types;all_psycho_struct.stim_vars]');
    b(1).FaceColor=[0 0 0];b(2).FaceColor=[.5 .5 .5];
    box off
    grid on
    ylabel('Indices')
    legend('Stim Type','Stim Var')
    xlabel('Trial types')
    xlim([0 num_trial_types+1])
end

%% get stim triggered average STAs
for i = 1:num_cells
    this_cell_trace = cnm_struct(cell_struct(i).cnm_idx).deconvC_full;
    if opt.IF_USE_FILTC
        this_online_trace = cell_struct(i).filtC; % online trace med filtered
    else
        this_online_trace = cnm_struct(cell_struct(i).cnm_idx).onlineC_full; % online trace without median filter
    end

    this_num_trials = numel(cnm_struct(cell_struct(i).cnm_idx).stim_frames );
    this_sens_stim_frames =  cnm_struct(cell_struct(i).cnm_idx).stim_frames;
    cell_struct(i).num_trials = this_num_trials*cell_struct(i).opsin_positive;
    cell_struct(i).is_sensory = 0;
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

% sta of backgound component (for alignment check)
[~,~,~,~,~,bg_sta_traces,bg_sta_trace] =...
    make_sta_from_traces(backgroundC,sens_stim_frames+opt.sta_stimon_frame ,opt.sta_pre_frames,opt.sta_post_frames,1:opt.sta_baseline_frames);
disp('got cell_struct sta')

%% Normlalise traces to baseline and std
% get baseline and std from the first sets of easy trials
easy_trial_idx = 1:10;
cell_bs = cell2mat(arrayfun(@(x)mean(reshape(cell_struct(x).raw_sta_traces(easy_trial_idx,1:opt.sta_baseline_frames),[1,numel(easy_trial_idx)*opt.sta_baseline_frames]),2),1:num_cells,'UniformOutput', false));
cell_sd = cell2mat(arrayfun(@(x)std(reshape(cell_struct(x).raw_sta_traces(easy_trial_idx,:),[1,numel(easy_trial_idx)*opt.trial_length])),1:num_cells,'UniformOutput', false));
for i = 1:num_cells
    
    cell_struct(i).raw_sta_traces = (cell_struct(i).raw_sta_traces - cell_bs(i))./cell_sd(i);
    cell_struct(i).raw_sta_trace = (cell_struct(i).raw_sta_trace - cell_bs(i))./cell_sd(i);
end
disp('normalised cell_struct raw_sta_traces')

%% sort STAs by different trial types
trial_types = fields(trial_indices);
raw_cell_struct = struct();
for i = 1:numel(trial_types)
    this_fd = trial_types{i};
    this_idx = trial_indices.(this_fd);
    this_idx = this_idx(this_idx<tot_num_trials);
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
%% check trial-history dependence
trial_hist = struct();
pre_num_trials = 2;
post_num_trials = 1;
for i = 1:numel(trial_types)
    this_fd = trial_types{i};
    this_idx = trial_indices.(this_fd);
    this_idx = this_idx(this_idx<tot_num_trials-post_num_trials&this_idx>pre_num_trials);
    [~,~,~,trial_hist.stim_type.(this_fd)] = make_sta_from_traces(trials.stim_type,this_idx,pre_num_trials,post_num_trials,1);
    [~,~,~,trial_hist.stim_var.(this_fd)] = make_sta_from_traces(trials.stim_var,this_idx,pre_num_trials,post_num_trials,1);
    [~,~,~,trial_hist.outcome.(this_fd)] = make_sta_from_traces(trials.correct,this_idx,pre_num_trials,post_num_trials,1);  
    trial_hist.reward_type.(this_fd) = trial_hist.stim_type.(this_fd).*trial_hist.outcome.(this_fd);  

end

%% GET CELL IDENTITY
% stimulus AUC calculated from correct trials with highest contrast and amp
stim_fds = {'stim_1_var_9_correct','stim_2_var_9_correct'}; 
[cell_struct] = get_cell_auc(cell_struct,stim_fds,'correct_stimAUC',opt);
correct_stimulusAUC_zscore = extractfield(cell_struct,'correct_stimAUC_zscore');
correct_texAUC = extractfield(cell_struct,'correct_stimAUC');
cell_idx_struct = struct();
cell_idx_struct.all = 1:num_cells;
cell_idx_struct.stim2 = find(correct_stimulusAUC_zscore>opt.N& correct_texAUC>0.55); % cells prefering texture1 in correct trials
cell_idx_struct.stim1 = find(correct_stimulusAUC_zscore<-opt.N& correct_texAUC<0.45); % cells prefering texture2 in correct trials
cell_idx_struct.stim = unique([cell_idx_struct.stim1,cell_idx_struct.stim2]);

for c = 1:num_cells
    cell_struct(c).is_tuned = 0;
    cell_struct(c).pref = 0;
    if ismember(c,cell_idx_struct.stim)
        cell_struct(c).is_tuned = 1;
        if ismember(c,cell_idx_struct.stim1)
            cell_struct(c).pref = 1;
        else
            cell_struct(c).pref = 2;
        end
    end
end
disp('got stim auc')
%% get choice auc
choice_fds = {'stim_1_var_1_correct','stim_1_var_1_incorrect'}; 
[cell_struct] = get_cell_auc(cell_struct,choice_fds,'choiceAUC',opt);
choiceAUC_zscore = extractfield(cell_struct,'choiceAUC_zscore');
choiceAUC = extractfield(cell_struct,'choiceAUC');
cell_idx_struct.port2 = find(choiceAUC_zscore>opt.N& choiceAUC>0.55); % cells prefering port2
cell_idx_struct.port1 = find(choiceAUC_zscore<-opt.N& choiceAUC<0.45); % cells prefering port1
cell_idx_struct.port = unique([cell_idx_struct.port1,cell_idx_struct.port2]);
disp('got choice auc')

%% get cells corelated with relative change - skip this section for online analysis
% tuned for the rewarded stim type
cmp_fds = {'stim_1_var_4_correct','stim_1_var_5_correct'}; % stim2 same, stim1 different
[cell_struct] = get_cell_auc(cell_struct,cmp_fds,'stim1AUC',opt);
cmp_fds = {'stim_2_var_2_correct','stim_2_var_4_correct'}; % stim1 same, stim2 different
[cell_struct] = get_cell_auc(cell_struct,cmp_fds,'stim2AUC',opt);
stim1AUC_zscore = extractfield(cell_struct,'stim1AUC_zscore');
stim2AUC_zscore = extractfield(cell_struct,'stim2AUC_zscore');
cell_idx_struct.stim1_tuned = find(abs(stim1AUC_zscore)>opt.N); 
cell_idx_struct.stim2_tuned = find(abs(stim2AUC_zscore)>opt.N);

% tuned for the unrewarded stim type
cmp_fds = {'stim_1_var_5_correct','stim_1_var_7_correct'}; % stim1 same, stim2 different
[cell_struct] = get_cell_auc(cell_struct,cmp_fds,'stim2AUC_oppo',opt);
cmp_fds = {'stim_2_var_4_correct','stim_2_var_5_correct'}; % stim2 same, stim1 different
[cell_struct] = get_cell_auc(cell_struct,cmp_fds,'stim1AUC_oppo',opt);
stim1AUC_oppo_zscore = extractfield(cell_struct,'stim1AUC_oppo_zscore');
stim2AUC_oppo_zscore = extractfield(cell_struct,'stim2AUC_oppo_zscore');
cell_idx_struct.stim1_oppo_tuned = find(abs(stim1AUC_oppo_zscore)>opt.N); 
cell_idx_struct.stim2_oppo_tuned = find(abs(stim2AUC_oppo_zscore)>opt.N);

%% discard cells with low snr or very few pixels
cell_snr = cell_bs./cell_sd.^2;
cell_size = cell2mat(arrayfun(@(x)numel(x.pix_values),cell_struct,'UniformOutput', false));

high_snr_cell_idx = find(cell_snr>opt.min_cell_snr&cell_size>opt.min_cell_size);
for  c = 1:num_cells
    cell_struct(c).snr = cell_snr(c);
    cell_struct(c).size = cell_size(c);
end
% delete low-snr cells from cell_idx_struct
fds = fields(cell_idx_struct);
for f = 1:numel(fds)
    fd = fds{f};
    cell_idx_struct.(fd) = intersect(cell_idx_struct.(fd),high_snr_cell_idx);
end

%% check high snr cells in FOV map
figure('name','cell snr check')
subplot(1,2,1)
imagesc(accepted_com_fov)
colormap(gray)
axis square
title('Accepted ROIs')
ax = subplot(1,2,2);
plot_value_in_rois( cell_struct, 'snr',[256 256],ax,'IF_NORM_PIX',0,...
    'IF_CONTOUR',0,'IF_SHOW_OPSIN',0,'target_cell_idx',cell_idx_struct.all);

%% MAKE OUTPUT FILE FOR PYRTAOI PHOTOEXCITABILITY TEST
[output1,save_time1] = generate_cell_idx_file(cell_struct,cell_idx_struct,[],opt);

%% Response amplitude vs ipsi deflection increase - skip this if not using a range of ipsi amplitudes
fds = {'stim_1_var_1_correct', 'stim_1_var_4_correct','stim_2_var_5_correct','stim_2_var_8_correct'};
lickright_fds = {'stim_1_var_1_incorrect', 'stim_1_var_4_incorrect','stim_2_var_5_correct','stim_2_var_8_correct'};

this_names = cellfun(@(x)strsplit(x,'_'),fds,'UniformOutput',false);
this_stim_types = cellfun(@(x)str2double(x(2)),this_names,'UniformOutput',false);
this_stim_vars = cellfun(@(x)str2double(x(4)),this_names,'UniformOutput',false);
[left,right]=cellfun(@(x,y) get_volts(x,y),this_stim_types,this_stim_vars);
for i = 1:num_cells
cell_struct(i).ipsi_correct_mean = cell2mat(cellfun(@(x)mean(mean(cell_struct(i).(x)(opt.sta_avg_frames,:),1)),fds,'UniformOutput',false));
cell_struct(i).ipsi_correct_std = cell2mat(cellfun(@(x)std(mean(cell_struct(i).(x)(opt.sta_avg_frames,:),1)),fds,'UniformOutput',false));
[this_corr,this_p] = corrcoef(cell_struct(i).ipsi_correct_mean',right');
cell_struct(i).ipsi_corr = this_corr(1,2);
cell_struct(i).ipsi_p = this_p(1,2);

cell_struct(i).ipsi_incorrect_mean= cell2mat(cellfun(@(x)mean(mean(cell_struct(i).(strrep(x,'_correct','_incorrect'))(opt.sta_avg_frames,:),1)),fds,'UniformOutput',false));
cell_struct(i).ipsi_incorrect_std = cell2mat(cellfun(@(x)std(mean(cell_struct(i).(strrep(x,'_correct','_incorrect'))(opt.sta_avg_frames,:),1)),fds,'UniformOutput',false));

cell_struct(i).ipsi_lickright_mean = cell2mat(cellfun(@(x)mean(mean(cell_struct(i).(x)(opt.sta_avg_frames,:),1)),lickright_fds,'UniformOutput',false));
cell_struct(i).ipsi_lickright_std = cell2mat(cellfun(@(x)mean(mean(cell_struct(i).(x)(opt.sta_avg_frames,:),1)),lickright_fds,'UniformOutput',false));
[this_corr,this_p] = corrcoef(cell_struct(i).ipsi_lickright_mean',right');
cell_struct(i).ipsi_lickright_corr = this_corr(1,2);
cell_struct(i).ipsi_lickright_p = this_p(1,2);

end
cell_idx_struct.ipsi_tuned =intersect( find(cell2mat({cell_struct(:).ipsi_p})<0.1),high_snr_cell_idx);

%% plot response amplitude aganst ipsi-deflection amplitude 
figure('name','cells tuned by ipsi input','units','normalized','outerposition',[0 0 1 1])
num_plot_cols = 8;
num_plot_rows = ceil(num_cells/num_plot_cols);
for i = 1:num_cells
    subtightplot(num_plot_rows,num_plot_cols,i)
    hold on
    errorbar(right,cell_struct(i).ipsi_correct_mean,cell_struct(i).ipsi_correct_std,'color','black')
    errorbar(right+0.1,cell_struct(i).ipsi_incorrect_mean,cell_struct(i).ipsi_incorrect_std,'color',[.7 .7 .7])

    xlim([-0.5 3])
    xlabel('Ipsi amp.')
    ylabel('Cell amp.')
    axis square
    title(['Cell' num2str(i),' r= ' num2str(cell_struct(i).ipsi_corr,'%10.2f')],'color',(cell_struct(i).ipsi_p<0.1).*[1 0 0])
end
export_fig([fig_save_path filesep 'IpsiTuning_' strrep(caiman_file,'.mat','.png')])

% plot ipsi tuning against AUCs
idx = cell_idx_struct.all;
corr_fd = 'ipsi_corr'; % high ipsi_corr means the cell is positively tuned by ipsi deflection, or overall deflection amplitude
choiceauc_fd = 'choiceAUC_zscore'; % high choice auc means the cell prefers lickport2(ipsi port)
stimauc_fd = 'correct_stimAUC_zscore';
figure('name','AUC and ipsi tuning')
subplot(2,2,1)
hold on
xx = cell2mat({cell_struct(idx).(stimauc_fd)});
yy = cell2mat({cell_struct(idx).(corr_fd)});
x = cell2mat({cell_struct(cell_idx_struct.ipsi_tuned ).(stimauc_fd)});
y = cell2mat({cell_struct(cell_idx_struct.ipsi_tuned ).(corr_fd)});
scatter(xx,yy,'MarkerEdgeColor','black')
scatter(x,y,'MarkerEdgeColor','red')
[ r,~,~,~,~,~,p ] = plot_fit_linear( xx,yy,[-5 5],'r');
title(['r=' num2str(r,'%10.2f') '; P= ' num2str(p,'%10.4f')])
axis square
xlabel('Stim AUC')
ylabel('Ipsi corr')

subplot(2,2,2)
hold on
xx = cell2mat({cell_struct(idx).(choiceauc_fd)});
yy = cell2mat({cell_struct(idx).(corr_fd)});
x = cell2mat({cell_struct(cell_idx_struct.ipsi_tuned ).(choiceauc_fd)});
y = cell2mat({cell_struct(cell_idx_struct.ipsi_tuned ).(corr_fd)});
scatter(xx,yy,'MarkerEdgeColor','black')
scatter(x,y,'MarkerEdgeColor','red')
[ r,~,~,~,~,~,p ] = plot_fit_linear( xx,yy,[-5 5],'r');
title(['r=' num2str(r,'%10.2f') '; P= ' num2str(p,'%10.4f')])
axis square
xlabel('Choice AUC')
ylabel('Ipsi corr')
subplot(2,2,3)
hold on
xx = abs(cell2mat({cell_struct(idx).(stimauc_fd)}));
yy = abs(cell2mat({cell_struct(idx).(corr_fd)}));
x = abs(cell2mat({cell_struct(cell_idx_struct.ipsi_tuned ).(stimauc_fd)}));
y = abs(cell2mat({cell_struct(cell_idx_struct.ipsi_tuned ).(corr_fd)}));
scatter(xx,yy,'MarkerEdgeColor','black')
scatter(x,y,'MarkerEdgeColor','red')
[ r,~,~,~,~,~,p ] = plot_fit_linear( xx,yy,[0 5],'r');
title(['r=' num2str(r,'%10.2f') '; P= ' num2str(p,'%10.4f')])

axis square
xlabel('Abs. Stim AUC')
ylabel('Abs. Ipsi corr')

subplot(2,2,4)
hold on
xx = abs(cell2mat({cell_struct(idx).(choiceauc_fd)}));
yy = abs(cell2mat({cell_struct(idx).(corr_fd)}));
x = abs(cell2mat({cell_struct(cell_idx_struct.ipsi_tuned ).(choiceauc_fd)}));
y = abs(cell2mat({cell_struct(cell_idx_struct.ipsi_tuned ).(corr_fd)}));
scatter(xx,yy,'MarkerEdgeColor','black')
scatter(x,y,'MarkerEdgeColor','red')
[ r,~,~,~,~,~,p ] = plot_fit_linear( xx,yy,[0 5],'r');
title(['r=' num2str(r,'%10.2f') '; P= ' num2str(p,'%10.4f')])
axis square
xlabel('Abs. Choice AUC')
ylabel('Abs. Ipsi corr')
export_fig([fig_save_path filesep 'IpsiTuningVsStimChoiceAUCs_' strrep(caiman_file,'.mat','.png')])


%% Compare AUCs - skip these for online analysis
figure('name','cells tuned by relative change','units','normalized','outerposition',[0 0 1 1])
subplot(2,4,1); hold on
histogram(stim1AUC_zscore,'FaceColor','none','EdgeColor',trial_color.L1,'DisplayStyle','stairs')
histogram(stim2AUC_zscore,'FaceColor','none','EdgeColor',trial_color.L2,'DisplayStyle','stairs')
legend('stim1 tuned','stim2 tuned')
axis square

subplot(2,4,2); hold on
scatter(stim1AUC_zscore,stim2AUC_zscore,'MarkerEdgeColor','black')
[ r,~,~,~,~,~,p ] = plot_fit_linear( stim1AUC_zscore,stim2AUC_zscore,[-5 5],'r');
axis square
xlabel('stim1 auc')
ylabel('stim2 auc')
title(['r=' num2str(r,'%10.2f') '; P= ' num2str(p,'%10.4f')])

subplot(2,4,3); hold on
scatter(correct_stimulusAUC_zscore,stim1AUC_zscore,'MarkerEdgeColor','black')
scatter(correct_stimulusAUC_zscore(cell_idx_struct.stim1_tuned),stim1AUC_zscore(cell_idx_struct.stim1_tuned),'MarkerEdgeColor',trial_color.L1)
scatter(correct_stimulusAUC_zscore(cell_idx_struct.stim2_tuned),stim1AUC_zscore(cell_idx_struct.stim2_tuned),'MarkerEdgeColor',trial_color.L2)
[ r,~,~,~,~,~,p ] = plot_fit_linear( correct_stimulusAUC_zscore,stim1AUC_zscore,[-5 5],'r');
axis square
xlabel('stim auc')
ylabel('stim1 auc')
title(['r=' num2str(r,'%10.2f') '; P= ' num2str(p,'%10.4f')])

subplot(2,4,4); hold on
scatter(correct_stimulusAUC_zscore,stim2AUC_zscore,'MarkerEdgeColor','black')
scatter(correct_stimulusAUC_zscore(cell_idx_struct.stim1_tuned),stim2AUC_zscore(cell_idx_struct.stim1_tuned),'MarkerEdgeColor',trial_color.L1)
scatter(correct_stimulusAUC_zscore(cell_idx_struct.stim2_tuned),stim2AUC_zscore(cell_idx_struct.stim2_tuned),'MarkerEdgeColor',trial_color.L2)

[ r,~,~,~,~,~,p ] = plot_fit_linear( correct_stimulusAUC_zscore,stim2AUC_zscore,[-5 5],'r');
axis square
xlabel('stim auc')
ylabel('stim2 auc')
title(['r=' num2str(r,'%10.2f') '; P= ' num2str(p,'%10.4f')])

subplot(2,4,5); hold on
histogram(stim1AUC_oppo_zscore,'FaceColor','none','EdgeColor',trial_color.L1,'DisplayStyle','stairs')
histogram(stim2AUC_oppo_zscore,'FaceColor','none','EdgeColor',trial_color.L2,'DisplayStyle','stairs')
legend('stim1 tuned','stim2 tuned')
axis square

subplot(2,4,6); hold on
scatter(stim1AUC_oppo_zscore,stim2AUC_oppo_zscore,'MarkerEdgeColor','black')
[ r,~,~,~,~,~,p ] = plot_fit_linear( stim1AUC_oppo_zscore,stim2AUC_oppo_zscore,[-5 5],'r');
axis square
xlabel('stim1 (oppo) auc')
ylabel('stim2 (oppo) auc')
title(['r=' num2str(r,'%10.2f') '; P= ' num2str(p,'%10.4f')])

subplot(2,4,7); hold on
scatter(correct_stimulusAUC_zscore,stim1AUC_oppo_zscore,'MarkerEdgeColor','black')
scatter(correct_stimulusAUC_zscore(cell_idx_struct.stim1_oppo_tuned),stim1AUC_oppo_zscore(cell_idx_struct.stim1_oppo_tuned),'MarkerEdgeColor',trial_color.L1)
scatter(correct_stimulusAUC_zscore(cell_idx_struct.stim2_oppo_tuned),stim1AUC_oppo_zscore(cell_idx_struct.stim2_oppo_tuned),'MarkerEdgeColor',trial_color.L2)
[ r,~,~,~,~,~,p ] = plot_fit_linear( correct_stimulusAUC_zscore,stim1AUC_oppo_zscore,[-5 5],'r');
axis square
xlabel('stim auc')
ylabel('stim1 (oppo) auc')
title(['r=' num2str(r,'%10.2f') '; P= ' num2str(p,'%10.4f')])

subplot(2,4,8); hold on
scatter(correct_stimulusAUC_zscore,stim2AUC_oppo_zscore,'MarkerEdgeColor','black')
scatter(correct_stimulusAUC_zscore(cell_idx_struct.stim1_oppo_tuned),stim2AUC_oppo_zscore(cell_idx_struct.stim1_oppo_tuned),'MarkerEdgeColor',trial_color.L1)
scatter(correct_stimulusAUC_zscore(cell_idx_struct.stim2_oppo_tuned),stim2AUC_oppo_zscore(cell_idx_struct.stim2_oppo_tuned),'MarkerEdgeColor',trial_color.L2)

[ r,~,~,~,~,~,p ] = plot_fit_linear( correct_stimulusAUC_zscore,stim2AUC_oppo_zscore,[-5 5],'r');
axis square
xlabel('stim auc')
ylabel('stim2 (oppo) auc')
title(['r=' num2str(r,'%10.2f') '; P= ' num2str(p,'%10.4f')])
%
figure
subplot(1,2,1); hold on
scatter(stim1AUC_zscore,stim1AUC_oppo_zscore,'MarkerEdgeColor','black')
[ r,~,~,~,~,~,p ] = plot_fit_linear( stim1AUC_zscore,stim1AUC_oppo_zscore,[-5 5],'r');
axis square
xlabel('stim1 auc')
ylabel('stim1 oppo auc')
title(['r=' num2str(r,'%10.2f') '; P= ' num2str(p,'%10.4f')])

subplot(1,2,2); hold on
scatter(stim2AUC_zscore,stim2AUC_oppo_zscore,'MarkerEdgeColor','black')
[ r,~,~,~,~,~,p ] = plot_fit_linear( stim2AUC_zscore,stim2AUC_oppo_zscore,[-5 5],'r');
axis square
xlabel('stim2 auc')
ylabel('stim2 oppo auc')
title(['r=' num2str(r,'%10.2f') '; P= ' num2str(p,'%10.4f')])

%%compare stim and choice auc
figure
subplot(1,5,1); hold on
scatter(choiceAUC_zscore,correct_stimulusAUC_zscore,'MarkerEdgeColor','black')
[ r,~,~,~,~,~,p ] = plot_fit_linear( abs(choiceAUC_zscore),abs(correct_stimulusAUC_zscore),[-5 5],'r');
axis square
xlabel('choice auc')
ylabel('stim auc')
title(['r=' num2str(r,'%10.2f') '; P= ' num2str(p,'%10.4f')])


subplot(1,5,2); hold on
scatter(choiceAUC_zscore,stim1AUC_zscore,'MarkerEdgeColor','black')
[ r,~,~,~,~,~,p ] = plot_fit_linear( choiceAUC_zscore,stim1AUC_zscore,[-5 5],'r');
axis square
xlabel('choice auc')
ylabel('stim1 auc')
title(['r=' num2str(r,'%10.2f') '; P= ' num2str(p,'%10.4f')])

subplot(1,5,3); hold on
scatter(choiceAUC_zscore,stim1AUC_oppo_zscore,'MarkerEdgeColor','black')
[ r,~,~,~,~,~,p ] = plot_fit_linear( choiceAUC_zscore,stim1AUC_oppo_zscore,[-5 5],'r');
axis square
xlabel('choice auc')
ylabel('stim1 oppo auc')
title(['r=' num2str(r,'%10.2f') '; P= ' num2str(p,'%10.4f')])


subplot(1,5,4); hold on
scatter(choiceAUC_zscore,stim2AUC_oppo_zscore,'MarkerEdgeColor','black')
[ r,~,~,~,~,~,p ] = plot_fit_linear( choiceAUC_zscore,stim2AUC_oppo_zscore,[-5 5],'r');
axis square
xlabel('choice auc')
ylabel('stim2 oppo auc')
title(['r=' num2str(r,'%10.2f') '; P= ' num2str(p,'%10.4f')])

subplot(1,5,5); hold on
scatter(choiceAUC_zscore,stim2AUC_zscore,'MarkerEdgeColor','black')
[ r,~,~,~,~,~,p ] = plot_fit_linear( choiceAUC_zscore,stim2AUC_zscore,[-5 5],'r');
axis square
xlabel('choice auc')
ylabel('stim2 auc')
title(['r=' num2str(r,'%10.2f') '; P= ' num2str(p,'%10.4f')])


%% Plot STAs for all components
figure('name','baseline session stim sta traces','units','normalized','outerposition',[0 0 1 1])
plot_num_cells = num_cells;
num_plot_cols = 8;
avg_frame_range = opt.sta_avg_frames;
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
    
    % mark tuned cells
    if( cell_struct(i).is_tuned)
        box on
        this_color = trial_color.(['L' num2str(cell_struct(i).pref)]);
        set(gca,'XColor',this_color,'YColor',this_color,'linewidth',3)
    end
    
    % mark time
    plot([1,1].*opt.sta_stimon_frame, ylim,'color','black','linestyle',':')
    if ~ opt.flag_use_peak
        plot([avg_frame_range(1),avg_frame_range(end)],[0,0],'color',[.5 .5 .5],'linewidth',2)
    else
        plot([peak_frame_range(1),peak_frame_range(end)],[0,0],'color',[.5 .5 .5],'linewidth',2)
    end
    
    % show auc
    text(0.05,.8,['stim auc '  num2str(cell_struct(i).correct_stimAUC_zscore,'%10.2f') ],'units','normalized', 'horizontalalignment','left', 'color','black')
    text(0.05,.7,['choice auc '  num2str(cell_struct(i).choiceAUC_zscore,'%10.2f') ],'units','normalized', 'horizontalalignment','left', 'color','black')
    
    % modify x ticks
    xaxisvalues = [0:45:opt.sta_pre_frames+opt.sta_post_frames];
    xticks(xaxisvalues)
    xticklabels(arrayfun(@(x){num2str(x)},(xaxisvalues-opt.sta_trialon_frame)./opt.frame_rate))
    
end
suptitle(strrep(strrep(caiman_file,'.mat',''),'_',' '))
export_fig([fig_save_path filesep 'StimSTATrace_' strrep(caiman_file,'.mat','.png')])

figure('name','baseline session choice sta traces','units','normalized','outerposition',[0 0 1 1])
plot_num_cells = num_cells;
num_plot_cols = 8;
avg_frame_range = opt.sta_avg_frames;
num_plot_rows = ceil(plot_num_cells/num_plot_cols);
for i = 1:plot_num_cells
    subtightplot(num_plot_rows,num_plot_cols,i)
    hold on
    % correct trials
    plot(cell_struct(i).(choice_fds{1}),'color',trial_color.L1,'linewidth',1)
    plot(cell_struct(i).(choice_fds{2}),'color',trial_color.L2,'linewidth',1)

    xlim([0 opt.sta_pre_frames+opt.sta_post_frames+1])
    ylim([-2 10])
    axis square
    
    text(1,1,['Cell ' num2str(i) ' (ROI ' num2str(cell_struct(i).cnm_idx) ')'],'units','normalized','color','black','Horizontalalignment','right','VerticalAlignment','top')

    if(~isempty(find(opsin_positive_idx==i)))
        text(1,1,['Cell ' num2str(i)],'units','normalized','color','r','Horizontalalignment','right','VerticalAlignment','top')
    end
    box off
    
    % mark tuned cells
    if( cell_struct(i).is_tuned)
        box on
        this_color = trial_color.(['L' num2str(cell_struct(i).pref)]);
        set(gca,'XColor',this_color,'YColor',this_color,'linewidth',3)
    end
    
    % mark time
    plot([1,1].*opt.sta_stimon_frame, ylim,'color','black','linestyle',':')
    if ~ opt.flag_use_peak
        plot([avg_frame_range(1),avg_frame_range(end)],[0,0],'color',[.5 .5 .5],'linewidth',2)
    else
        plot([peak_frame_range(1),peak_frame_range(end)],[0,0],'color',[.5 .5 .5],'linewidth',2)
    end
    
    % show auc
    text(0.05,.8,['stim auc '  num2str(cell_struct(i).correct_stimAUC_zscore,'%10.2f') ],'units','normalized', 'horizontalalignment','left', 'color','black')
    text(0.05,.7,['choice auc '  num2str(cell_struct(i).choiceAUC_zscore,'%10.2f') ],'units','normalized', 'horizontalalignment','left', 'color','black')
    
    % modify x ticks
    xaxisvalues = [0:45:opt.sta_pre_frames+opt.sta_post_frames];
    xticks(xaxisvalues)
    xticklabels(arrayfun(@(x){num2str(x)},(xaxisvalues-opt.sta_trialon_frame)./opt.frame_rate))
    
end
suptitle(strrep(strrep(caiman_file,'.mat',''),'_',' '))
export_fig([fig_save_path filesep 'ChoiceStimSTATrace_' strrep(caiman_file,'.mat','.png')])


%% Train decoders
% add some fields to trialcolor
for i = 1:numel(trial_types)
    this_fd = trial_types{i};
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
opt.trial_color = trial_color;
%% manaully select threshold conditions - make it automatic later
TS_L1_fds = {'stim_1_var_1_correct'};
TS_L2_fds = {'stim_1_var_1_incorrect'};
opt.trigger_idx_fd = 'all';
%% factor analysis
% using trials with highest contrast plus the selected threshld fds
high_contrast_fds = {'stim_1_var_9_correct','stim_1_var_9_incorrect','stim_2_var_9_correct','stim_2_var_9_incorrect',...
                   'stim_2_var_6_correct','stim_2_var_6_incorrect',...
                   'stim_1_var_3_correct','stim_1_var_3_incorrect'};
fa_opt.bin_size = 1;
fa_opt.gocue_bin = floor(opt.sta_gocue_frame/fa_opt.bin_size);
fa_opt.stim_bin = ceil(opt.sta_stimon_frame/fa_opt.bin_size);
fa_opt.frame_rate = 30/fa_opt.bin_size;
fa_opt.Fs = opt.frame_rate;
fa_opt.trial_length = opt.trial_length/fa_opt.bin_size;
fa_opt.trial_color = trial_color;

fa_opt.idx_fields = {opt.trigger_idx_fd};
fa_opt.fd_names = unique([high_contrast_fds, TS_L1_fds,TS_L2_fds]);
fa_opt.plot_fds = fa_opt.fd_names;
fa_opt.m = 3;
fa_opt.IF_MEDFILT = 0;
disp('Rnning factor analysis...')
tic
[traces_in,fa_trial_idx,num_trials] = get_input_seq(raw_cell_struct,cell_idx_struct.(fa_opt.idx_fields{1}),...
    fa_opt.fd_names,fa_opt.bin_size,'IF_MEDFILT',fa_opt.IF_MEDFILT);%%
fa_struct = struct();
fa_struct.mean = mean(traces_in);
fa_struct.std = std(traces_in);
[fa_struct.lambda,fa_struct.psi,fa_struct.T,fa_struct.stats,fa_struct.F] = factoran(traces_in,fa_opt.m,'Xtype','data','Maxit',1000);
invsqrtPsi = diag(1 ./  sqrt(fa_struct.psi)); % get transition matrix (multiply to zscored data)
fa_struct.transmat = invsqrtPsi/(fa_struct.lambda'*invsqrtPsi);


fa_trial_idx.TS_L1 = cell2mat(cellfun(@(x)fa_trial_idx.(x),TS_L1_fds,'UniformOutput' , false));
fa_trial_idx.TS_L2 = cell2mat(cellfun(@(x)fa_trial_idx.(x),TS_L2_fds,'UniformOutput' , false));

[fa_traj_struct,num_fs] = get_pop_vectors(fa_struct.F,fa_opt.trial_length,fa_trial_idx);
toc
disp('...Done')
plot_pop_vectors(fa_traj_struct,fa_opt.plot_fds,fa_opt.m,fa_opt,...
    'plot_ylabel','Factor level','plot_num_cols',2);
%% hidden markov 
    hmm_opt = fa_opt;
    hmm_opt.plot_fds = hmm_opt.fd_names;
    hmm_opt.num_states = 3;
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
    
    disp('Rnning HMM analysis...')
    tic
    T = hmm_opt.trial_length * ones(num_trials,1);
    [hmm,Gamma_hmm, ~, vpath] = hmmmar(fa_struct.F,T,hmm_opt);
    toc
    disp('...Done')
    [hmm_struct,num_states] = get_pop_vectors(Gamma_hmm,hmm_opt.trial_length - hmm_opt.order,fa_trial_idx);
    hmm_opt.m = num_states;
    %plot states and paths
    plot_pop_vectors(hmm_struct,hmm_opt.plot_fds,num_states,hmm_opt,'IF_MEDIAN',0,'plot_num_cols',4, 'plot_area',true);
    vpath_struct = get_pop_vectors(vpath,hmm_opt.trial_length-hmm_opt.order,fa_trial_idx);
    plot_pop_vectors(vpath_struct,hmm_opt.plot_fds,1,hmm_opt,...
        'ylimit',[0, hmm_opt.K+1],'plot_ylabel','Viterbi path','IF_MEDIAN',0,'plot_num_cols',4);

%% get coding direction
% port cells in high contrast correct trials
% using raw onine traces to mimic online condition
coding_direct = cell2mat(arrayfun(@(x)mean(mean(cell_struct(x).(choice_fds{1})(opt.sta_avg_frames,:),2)),cell_idx_struct.(opt.trigger_idx_fd),'UniformOutput', false))...
    - cell2mat(arrayfun(@(x)mean(mean(cell_struct(x).(choice_fds{2})(opt.sta_avg_frames,:),2)),cell_idx_struct.(opt.trigger_idx_fd),'UniformOutput', false));
% project to coding direction
cd_traj_struct = get_cd_projection(cell_struct(cell_idx_struct.(opt.trigger_idx_fd)),coding_direct);
plot_pop_vectors(cd_traj_struct,choice_fds,1,opt,...
        'plot_ylabel','Projection')

%% choose which to use - factor analysis, or coding direction
if opt.IF_USE_FA
    traj_struct = fa_traj_struct;
elseif opt.IF_USE_HMM
    traj_struct = vpath_struct;
else
    traj_struct = cd_traj_struct;
end
%% get stim decoder (using easy trials)
stim_opt = fa_opt;
stim_opt.fd_names = stim_fds; % stim 1 will be positive
stim_opt.frames_to_avg = stim_opt.stim_bin+15:1:stim_opt.gocue_bin;
stim_opt.frames_to_train = round([1:1:150]/stim_opt.bin_size);
stim_opt.Nstd = 2;
stim_opt.min_frames = 15;
stim_opt.IF_FRAMEWISE =0;
stim_struct = {};
stim_proj_struct = {};
disp('Rnning stim decoder...')
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

% compare projection on stim axis trials
figure
values = struct();
values.stim1 = mean(stim_proj_struct.(stim_fds{1})(:,opt.sta_avg_frames),2) ;
values.stim2 = mean(stim_proj_struct.(stim_fds{2})(:,opt.sta_avg_frames),2) ;
scatter_cmp_conditions(values,[],...
    1,[0 0 0;0 0 0],'connect_scatter',0,'BriefXlabel',0,'ShowMeanInXlabel',1,'add_jitter',1);
axis square

%% get choice decoder (directly from traj, using threshold trials) 
choice_opt = fa_opt;
choice_opt.fd_names = {'TS_L2','TS_L1'}; % Choice 2 will be positive
choice_opt.fd_names = {TS_L2_fds{1},TS_L1_fds{1}};
choice_opt.frames_to_avg = [100:130]; %using baseline to decode choice
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
opt.pop_opt = stim_opt;
stim_thresh = -stim_struct.B(:,1);

if opt.IF_USE_FA
    stim_weights = fa_struct.transmat* stim_struct.B(:,2:end)';
    [stim_norm_weights,stim_norm_thresh] = get_norm_weights(stim_weights,stim_thresh,fa_struct.mean,fa_struct.std); % called in CompareDecoderWithAnimal
%     stim_norm_weights = stim_weights;
%     stim_norm_thresh = stim_thresh;

else
    stim_weights =coding_direct'* stim_struct.B(2:end)';
    stim_norm_weights = stim_weights;
    stim_norm_thresh = stim_thresh;
end

test_fd_names = {'stim_1_var_3_correct','stim_1_var_3_incorrect','stim_1_var_3_miss',...
                     'stim_1_var_4_correct','stim_1_var_4_incorrect','stim_1_var_4_miss',...
                     'stim_1_var_5_correct','stim_1_var_5_incorrect','stim_1_var_5_miss',...
                     'stim_1_var_7_correct','stim_1_var_7_incorrect','stim_1_var_7_miss',...
                     'stim_2_var_2_correct','stim_2_var_2_incorrect','stim_2_var_2_miss'...
                     'stim_2_var_4_correct','stim_2_var_4_incorrect','stim_2_var_4_miss'...
                     'stim_2_var_5_correct','stim_2_var_5_incorrect','stim_2_var_5_miss'...
                     'stim_2_var_6_correct','stim_2_var_6_incorrect','stim_2_var_6_miss'};

stim_proj_struct = struct();
[stim_proj_struct] = get_projections(raw_cell_struct(cell_idx_struct.(opt.trigger_idx_fd)),stim_norm_weights,test_fd_names,'proj_struct',stim_proj_struct,'bias',-stim_norm_thresh,'IS_CELL_STRUCT',1);

plot_pop_vectors(stim_proj_struct,test_fd_names,1,opt,...
        'ylimit',[-10 10],'plot_ylabel','Projection','plot_num_cols',6,'IF_PLOT_RAW_ONLY',1)
suptitle('Stim decoder projections')
% export_fig([fig_save_path filesep 'StimDecodProject_' strrep(caiman_file,'.mat','.png')])
%% get projections on choice decoder
choice_thresh = -choice_struct.B(:,1);
if opt.IF_USE_FA
    choice_weights = fa_struct.transmat* choice_struct.B(:,2:end)';
    [choice_norm_weights,choice_norm_thresh] = get_norm_weights(choice_weights,choice_thresh,fa_struct.mean,fa_struct.std); % called in CompareDecoderWithAnimal

else
    choice_norm_weights = coding_direct'* choice_struct.B(:,2:end)';
    choice_norm_thresh = choice_thresh;
end

choice_proj_struct = struct();
[choice_proj_struct] = get_projections(raw_cell_struct(cell_idx_struct.(opt.trigger_idx_fd)),choice_norm_weights,...
    test_fd_names,'proj_struct',choice_proj_struct,'bias',-choice_norm_thresh,'IS_CELL_STRUCT',1);

plot_pop_vectors(choice_proj_struct,test_fd_names,1,opt,...
    'ylimit',[-1000 1000],'plot_ylabel','Projection','plot_num_cols',6,'IF_PLOT_RAW_ONLY',1)
suptitle('Choice decoder projections')


%% get choice decoder after projecting to the stim axis
% L2 comes first
% keep IF_CROSSVAL = 1 if sample size is not too small(<5), to avoid biased
% by single trials in the unshuffled field
decod_fds = { 'stim_1_var_5_incorrect','stim_1_var_5_correct'}; % first fd will be positive; pyRTAOI takes positive stim1 value as incorrect trial
choice_after_stim_struct = struct();
choice_after_stim_struct =  get_binary_classifier( choice_after_stim_struct,stim_proj_struct, choice_opt,...
    'IF_CROSSVAL',1,'IF_FRAMEWISE',0,'fd_names',decod_fds);
choicestim_proj_struct = struct();
 [choicestim_proj_struct] = get_projections(stim_proj_struct,choice_after_stim_struct.B(:,2:end)',decod_fds,'proj_struct',choicestim_proj_struct,'bias',choice_after_stim_struct.B(:,1));
    
plot_pop_vectors(choicestim_proj_struct,decod_fds,1,opt,...
    'ylimit',[-50 50],'plot_ylabel','Projection','plot_num_cols',2,'IF_PLOT_RAW_ONLY',0)

%% get weights and thresh after the two decoders
norm_weights = stim_norm_weights*choice_after_stim_struct.B(2);
norm_thresh = -((-stim_norm_thresh*choice_after_stim_struct.B(2))+choice_after_stim_struct.B(1));
weights = stim_weights*choice_after_stim_struct.B(2);
thresh = -((-stim_thresh*choice_after_stim_struct.B(2))+choice_after_stim_struct.B(1));
% weights = fa_struct.transmat* stim_struct.B(2:end)'*choice_after_stim_struct.B(2);
% thresh = -((-stim_struct.thresh_fix*choice_after_stim_struct.B(2))+choice_after_stim_struct.B(1));

[choicestim_proj_struct] = get_projections(raw_cell_struct(cell_idx_struct.(opt.trigger_idx_fd)),norm_weights,test_fd_names,'bias',-norm_thresh,'IS_CELL_STRUCT',1);

[ choice_after_stim_struct ] =  get_binary_decoder_disc_time(choicestim_proj_struct, choice_after_stim_struct, ...
    decod_fds,choice_opt,'IF_FRAMEWISE',choice_opt.IF_FRAMEWISE,'threshold',0,'IF_USE_CROSSVAL_RESULT',false);
plot_pop_vectors(choicestim_proj_struct,decod_fds,1,opt,...
        'ylimit',[-1 1],'plot_ylabel','Projection','IF_PLOT_RAW_ONLY',0)

figure; 
hold on
plot_binary_decoder(choice_after_stim_struct,choice_opt)
suptitle('Choice-stim decoder')
% export_fig([fig_save_path filesep 'ChoiceStimDecoderPerform_' strrep(caiman_file,'.mat','.png')])

%% get projections on choice-stim decoder
plot_pop_vectors(choicestim_proj_struct,test_fd_names,1,opt,...
       'ylimit',[-.5 .5],'plot_ylabel','proj to choice-stim','plot_num_cols',3,'IF_PLOT_RAW_ONLY',1)
suptitle('Choice-stim decoder')
% export_fig([fig_save_path filesep 'ChoiceStimDecodProject_' strrep(caiman_file,'.mat','.png')])

%% plot accuracy - select tiral types to condition according to these plots
test_decod_struct = [];
num_compares = round(numel(test_fd_names)/3);
for i = 1:num_compares
    IF_REVERSE = 1;
    this_correct_fd = test_fd_names{3*(i-1)+1};
    this_incorrect_fd = test_fd_names{3*(i-1)+2};
    this_cmp_fds = {this_correct_fd,this_incorrect_fd};
    this_stim_type = strsplit(this_correct_fd,'_');
    this_stim_type = this_stim_type{2};
    % reverse for stim type2
    if strcmp(this_stim_type,'2')||strcmp(this_stim_type,'3')
        IF_REVERSE = 0;
    end
    try
        [  test_decod_struct{i} ] =  get_binary_decoder_disc_time( choicestim_proj_struct, choice_after_stim_struct,...
            this_cmp_fds,choice_opt,'IF_FRAMEWISE',choice_opt.IF_FRAMEWISE,'threshold',0,'IF_REVERSE',IF_REVERSE);
        test_decod_struct{i}.correct_fd = this_correct_fd;
        test_decod_struct{i}.incorrect_fd = this_incorrect_fd;
        
        figure('name','decoder performance tests')
        plot_binary_decoder(test_decod_struct{i},choice_opt)
        title(['No. trials: ' num2str(size(choicestim_proj_struct.(this_correct_fd),1)) ', ' num2str(size(choicestim_proj_struct.(this_incorrect_fd),1))])
        suptitle(cellfun(@(x)strrep(x,'_',' '),this_cmp_fds, 'UniformOutput', false))
        
    catch
        test_decod_struct{i} = [];
    end
end
%% plot accuracy of choice decoder
test_decod_struct = [];
num_compares = round(numel(test_fd_names)/3);
for i = 1:num_compares
    IF_REVERSE = 1;
    this_correct_fd = test_fd_names{3*(i-1)+1};
    this_incorrect_fd = test_fd_names{3*(i-1)+2};
    this_cmp_fds = {this_correct_fd,this_incorrect_fd};
    this_stim_type = strsplit(this_correct_fd,'_');
    this_stim_type = this_stim_type{2};
    % reverse for stim type2
    if strcmp(this_stim_type,'2')||strcmp(this_stim_type,'3')
        IF_REVERSE = 0;
    end
    try
        [  test_decod_struct{i} ] =  get_binary_decoder_disc_time( choice_proj_struct, choice_struct,...
            this_cmp_fds,choice_opt,'IF_FRAMEWISE',choice_opt.IF_FRAMEWISE,'threshold',0,'IF_REVERSE',IF_REVERSE);
        test_decod_struct{i}.correct_fd = this_correct_fd;
        test_decod_struct{i}.incorrect_fd = this_incorrect_fd;
        
        figure('name','decoder performance tests')
        plot_binary_decoder(test_decod_struct{i},choice_opt)
        title(['No. trials: ' num2str(size(choice_proj_struct.(this_correct_fd),1)) ', ' num2str(size(choice_proj_struct.(this_incorrect_fd),1))])
        suptitle(cellfun(@(x)strrep(x,'_',' '),this_cmp_fds, 'UniformOutput', false))
        
    catch
        test_decod_struct{i} = [];
    end
end
%% SELECT CONDITION STIM TYPES
disp('ENTER HERE!')
decod_struct = choice_struct;
condition_stim_type = [ 1 1];
condition_stim_var  = [ 1 1];
condition_type = {[101],[101]};
opt.target_idx_fd = {'all','all'};
pop_weights = choice_norm_weights; % select which weights to use
pop_thresh = choice_norm_thresh;

% generate parameters for pyRTAOI population analysis
pop_params = struct();
pop_params.weights = pop_weights;
pop_params.thresh = pop_thresh;
pop_params.frames_enable_trigger = max(opt.sta_stimon_frame-opt.sta_baseline_frames+15,fa_opt.bin_size*(decod_struct.shuf_disc_frame-opt.sta_baseline_frames))+[0 15];
pop_params.condition_stim_type = condition_stim_type; % CHANGE THIS
pop_params.condition_stim_var  = condition_stim_var; % type and var need to match
pop_params.condition_type = condition_type; % match target ensembles 100*stim_type + stim_var
frames_of_interest = [pop_params.frames_enable_trigger(1):pop_params.frames_enable_trigger(end)]+opt.sta_baseline_frames;
% [proj_sd,fds_of_interest] = get_proj_noise(choice_proj_struct,condition_stim_type,condition_stim_var,frames_of_interest);
% plot_pop_vectors(choice_proj_struct,fds_of_interest,1,opt,...
%        'noise_thresh',proj_sd,'ylimit',[-100 100],'plot_ylabel','proj to choice-stim','plot_num_cols',2,'IF_PLOT_RAW_ONLY',1)
pop_params.thresh_sd = 0;% NOT USING SD
% pop_params.fds_of_interest = fds_of_interest;
pop_params.frames_enable_trigger  = [110 125]-30;
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
    cell_idx_struct.photo_stim1 = intersect(photo_idx.all, cell_idx_struct.stim1);
    cell_idx_struct.photo_stim2 = intersect(photo_idx.all, cell_idx_struct.stim2);
    opt.target_idx_fd = {'photo_stim1','photo_stim2'}; % overide target indices with the excitable ones
end
% opt.target_idx_fd = {'stim1','stim2'};

%% =================== MAKE OUTPUT FILE FOR PYRTAOI =======================
[output2,save_time2] = generate_cell_idx_file(cell_struct,cell_idx_struct,pop_params,opt);

%% ==================== MAKE PYBEHAV EXTERNAL FILE ========================
condition_pybehav_stimtypes = [5,4]; % these stim types will not be rewarded
loop_stim_types = [arrayfun(@(x)condition_pybehav_stimtypes(x),pop_params.condition_stim_type) 1 2]; % condition trial types plus two extreme stim types
loop_stim_vars  = [pop_params.condition_stim_var  3 6];
pyrtaoi_stim_types = [pop_params.condition_stim_type 1 2]; 
num_loops = 30;
file_save_name = [opt.exp_name '_PyBehavior_' save_time2];
[trial_seq] = generate_deflection_trial_seq(loop_stim_types,loop_stim_vars,num_loops,opt.output_path,file_save_name);
temp_types = trial_seq(1,:);
temp_types(temp_types==condition_pybehav_stimtypes(1))=1;
temp_types(temp_types==condition_pybehav_stimtypes(2))=2;
pyrtaoi_seq = trial_seq;
pyrtaoi_seq(1,:) = temp_types;
pyrtaoi_file_save_name = [opt.exp_name '_RTAOiPyBehavior_' save_time2];
dlmwrite([opt.output_path filesep pyrtaoi_file_save_name '.txt'],pyrtaoi_seq)
disp(['saved as:' opt.output_path filesep pyrtaoi_file_save_name '.txt'] )

%% save cell structure 
output_save_name = [save_path filesep  'ProcDeflect_' caiman_file];
save(output_save_name,'cell_struct')
disp(['Output struct saved as:' output_save_name ])


%% ============================     END    ================================




%% stuff tried but didnt help
%% loop through all stimulus types and try to decode choice using projection onto stimulus direction
% - didnt help
all_decod_struct ={};
all_decod_proj = struct();
num_compares = round(numel(test_fd_names)/3);
for i = 1:num_compares
    IF_REVERSE = 1;
    this_correct_fd = test_fd_names{3*(i-1)+1};
    this_incorrect_fd = test_fd_names{3*(i-1)+2};
    this_cmp_fds = {this_correct_fd,this_incorrect_fd};
    this_stim_type = strsplit(this_correct_fd,'_');
    this_stim_type = this_stim_type{2};
    % reverse for stim type2
    if strcmp(this_stim_type,'2')||strcmp(this_stim_type,'3')
        IF_REVERSE = 0;
    end
    try
        all_decod_struct{i} = get_binary_classifier( struct(),stim_proj_struct, choice_opt,...
            'IF_CROSSVAL',0,'IF_FRAMEWISE',0,'fd_names',this_cmp_fds);
        [all_decod_proj] = get_projections(stim_proj_struct,all_decod_struct{i}.B(:,2:end)',this_cmp_fds,'proj_struct',all_decod_proj,'bias', all_decod_struct{i} .B(:,1));

        [  all_decod_struct{i} ] =  get_binary_decoder_disc_time( all_decod_proj, all_decod_struct{i},...
            this_cmp_fds,choice_opt,'IF_FRAMEWISE',choice_opt.IF_FRAMEWISE,'threshold',0,'IF_REVERSE',IF_REVERSE);
        all_decod_struct{i}.correct_fd = this_correct_fd;
        all_decod_struct{i}.incorrect_fd = this_incorrect_fd;
        
        figure('name','decoder performance tests')
        plot_binary_decoder(all_decod_struct{i},choice_opt)
        title(['No. trials: ' num2str(size(choice_proj_struct.(this_correct_fd),1)) ', ' num2str(size(choice_proj_struct.(this_incorrect_fd),1))])
        suptitle(cellfun(@(x)strrep(x,'_',' '),this_cmp_fds, 'UniformOutput', false))
        
    catch
        all_decod_struct{i} = [];
    end
end
