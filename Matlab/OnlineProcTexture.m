% generate sta from pyRTAOI result file after a sensory stim movie
% note: rois idx doesn't match with what's shown in pyRTAOI image window
% run this to quickly find sensory-opsin-positve cells  

% ZZ 2018

% adapted from OnlineProcVisual

% TO DO:
% concatinate multiple sessions
% check different epochs; classify cells baseed on epoch time-average
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
% matlab_set_paths_zz

%% parameters - CHANGE THIS
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

opt.correct_trial_only = true; % only use correct trials to get tunning

% select which online trace to use
opt.IF_USE_FILTC = true;

% select which dim reduction to use
opt.IF_USE_FA = true;

[trial_color] = online_tex_init_color();
%% load CAIMAN data
[caiman_file,caiman_path] = uigetfile('*.mat','Select texture caiman data');
disp(['Loaded file :',fullfile(caiman_path,caiman_file)])
caiman_data = load(fullfile(caiman_path,caiman_file)); 
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
%% organise data (generate plots for sanity check)
tex_stim_frames = {};
if(~isempty(caiman_data.stim_frames_caiman))
    sens_stim_frames = caiman_data.sensory_stim_frames+caiman_data.t_init;
else
    sens_stim_frames = [];
end

% discard sens stim frames beyond tottime (deal with senarios when movie is aborted early)
discard_trial_idx = find(sens_stim_frames>caiman_data.t_cnm-opt.sta_post_frames-opt.sta_pre_frames);
sens_stim_frames(discard_trial_idx) = [];
num_trials = length(sens_stim_frames);
fds = fields(trials);
for f =1:numel(fds)
trials.(fds{f})(discard_trial_idx) = [];
end
% trial type
trialOrder = caiman_data.trialOrder;
trialOrder(discard_trial_idx) = [];
trialTypes = unique(trialOrder);
num_stim_type = length(unique(trialOrder)); % orientations/texture

% only get correct trials
if opt.correct_trial_only &&  FLAG_PYBEHAV_LOADED
%     sens_stim_frames = sens_stim_frames(trials.correct==1);
%     trialOrder = trialOrder(trials.correct==1);
else
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

num_trials = numel(sens_stim_frames);
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
    temp_trace(caiman_frames) =  caiman_data.noisyC(i+1,1:num_frames);
    temp_trace = fillmissing(temp_trace,'linear');
    cnm_struct(i).noisyC_full = temp_trace;
    
    % use go-cue frame as 'stim frame' for sta traces
    cnm_struct(i).stim_frames = sens_stim_frames+opt.sta_gocue_frame;

end

temp_trace = nan(1,tot_frames);
temp_trace(caiman_frames) =  caiman_data.noisyC(1,1:num_frames);
temp_trace = fillmissing(temp_trace,'linear');
backgroundC = temp_trace;

%% only get accepted cells
accepted_idx = caiman_data.accepted_idx+1;
num_cells = numel(accepted_idx);

%% only analyse frames from the current recording
glob_trialon_frames = caiman_data.sensory_stim_frames + caiman_data.t_init;
% glob_trialon_frames(glob_trialon_frames>caiman_data.t_cnm-opt.sta_post_frames-opt.sta_pre_frames) = [];
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
    
    cell_struct(i).filtC = caiman_data.filt_C(i,1:num_frames);

    if(~isempty(find(opsin_positive_idx==i)))
         cell_struct(i).opsin_positive = 1;
    end
end

%% plot traces
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
yticks([ 1:num_comp].*plot_offset)
yticklabels(1:num_comp)
xlim([caiman_data.t_init tot_frames])
ylim([0 num_comp].*plot_offset+5)

% background
plot(backgroundC,'color',[.5 .5 .5],'linestyle',':')

for i = 1:numel(glob_trialon_frames)
    this_color = trial_color.(['stim' num2str(trials.stim_type(i) )]);
    plot([glob_trialon_frames(i) glob_trialon_frames(i)],ylim,'color',this_color) % trial-on 
    plot([glob_trialon_frames(i) glob_trialon_frames(i)]+opt.gocue_frame,ylim,'color',this_color,'linestyle',':') % go-cue
end

%% stim triggered average 
for i = 1:num_cells
    this_cell_trace = cnm_struct(cell_struct(i).cnm_idx).deconvC_full;
    
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
% this is implemented for online 
easy_trial_idx = 1:10;
cell_bs = cell2mat(arrayfun(@(x)mean(reshape(cell_struct(x).raw_sta_traces(easy_trial_idx,1:opt.sta_baseline_frames),[1,numel(easy_trial_idx)*opt.sta_baseline_frames]),2),1:num_cells,'UniformOutput', false));
cell_sd = cell2mat(arrayfun(@(x)std(reshape(cell_struct(x).raw_sta_traces(easy_trial_idx,:),[1,numel(easy_trial_idx)*opt.trial_length])),1:num_cells,'UniformOutput', false));
for i = 1:num_cells
    
    cell_struct(i).raw_sta_traces = (cell_struct(i).raw_sta_traces - cell_bs(i))./cell_sd(i);
    cell_struct(i).raw_sta_trace = (cell_struct(i).raw_sta_trace - cell_bs(i))./cell_sd(i);
end
disp('normalised cell_struct raw_sta_traces')


%% Get trial types
trial_indices = struct(); % % get trialtype-outcome indices
for i = 1:num_stim_type
    trial_indices.(['stim_' num2str(i) '_correct']) = find(trials.correct==1&trials.stim_type==i&trials.cheated==0);
    trial_indices.(['stim_' num2str(i) '_incorrect']) = find(trials.incorrect==1&trials.stim_type==i&trials.cheated==0);
end

%% STAs for different trial types
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

%% GET CELL IDENTITY
% stimulus AUC calculated from correct trials
peak_frame_range = opt.sta_peak_search_range;
avg_frame_range = opt.sta_avg_frames;
[cell_struct] = get_cell_auc(cell_struct,{'stim_1_correct','stim_2_correct'},'correct_stimAUC',opt);
[cell_struct] = get_cell_auc(cell_struct,{'stim_1_incorrect','stim_2_incorrect'},'incorrect_stimAUC',opt);
for c = 1:num_cells
    cell_struct(c).is_tuned = 0;
    cell_struct(c).pref = 0;
    if ismember(c,cell_idx_struct.tex)
        cell_struct(c).is_tuned = 1;
        if ismember(c,cell_idx_struct.tex1)
            cell_struct(c).pref = 1;
        else
            cell_struct(c).pref = 2;
        end
    end
end
correct_stimulusAUC_zscore = extractfield(cell_struct,'correct_stimAUC_zscore');
correct_texAUC = extractfield(cell_struct,'correct_stimAUC');
%%
cell_idx_struct = struct();
cell_idx_struct.all = 1:num_cells;
cell_idx_struct.tex2 = find(correct_stimulusAUC_zscore>opt.N& correct_texAUC>0.55); % cells prefering texture1 in correct trials
cell_idx_struct.tex1 = find(correct_stimulusAUC_zscore<-opt.N& correct_texAUC<0.45); % cells prefering texture2 in correct trials
cell_idx_struct.tex = unique([cell_idx_struct.tex1,cell_idx_struct.tex2]);

incorrect_stimulusAUC_zscore = extractfield(cell_struct,'incorrect_stimAUC_zscore');
incorrect_texAUC = extractfield(cell_struct,'incorrect_stimAUC');
cell_idx_struct.port2 = find(correct_stimulusAUC_zscore>opt.N & incorrect_stimulusAUC_zscore<-opt.N); % cells prefering texture1 in correct trials
cell_idx_struct.port1 = find(correct_stimulusAUC_zscore<-opt.N& incorrect_stimulusAUC_zscore>opt.N); % cells prefering texture2 in correct trials
cell_idx_struct.port = unique([cell_idx_struct.port1,cell_idx_struct.port2]);

cell_idx_struct.relevant = unique([cell_idx_struct.port,cell_idx_struct.tex]);

disp('got cell identity')


%% discard cells with low snr or very few pixels (optional)
cell_snr = cell_bs./cell_sd.^2;
cell_size = cell2mat(arrayfun(@(x)numel(x.pix_values),cell_struct,'UniformOutput', false));
for  c = 1:num_cells
    cell_struct(c).snr = cell_snr(c);
    cell_struct(c).size = cell_size(c);
end
figure('name','cell snr check')
subplot(1,2,1)
imagesc(accepted_com_fov)
colormap(gray)
axis square
title('Accepted ROIs')
ax = subplot(1,2,2);
plot_value_in_rois( cell_struct, 'snr',[256 256],ax,'IF_NORM_PIX',0,...
    'IF_CONTOUR',0,'IF_SHOW_OPSIN',0,'target_cell_idx',cell_idx_struct.all);
title('SNR')
%% select cell identity for readout and stimulation
opt.target_idx_fd = {'tex1','tex2'};
opt.trigger_idx_fd = 'relevant';
opt.fov_size = double(cnm_dims);
opt.ds_factor = caiman_data.ds_factor;

%% MAKE OUTPUT FILE FOR PYRTAOI PHOTOEXCITABILITY TEST
[output1,save_time1] = generate_cell_idx_file(cell_struct,cell_idx_struct,[],opt);

%% Plot STAs traces for all components (save figures)
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

figure('name','condition sta traces (errobar)','units','normalized','outerposition',[0 0 1 1])
x_ticks =[0:1:opt.trial_length-1];

for i = 1:num_cells
    subtightplot(num_plot_rows,num_plot_cols,i)
    hold on
    for j = 1:num_stim_type
        this_traces = cell_struct(i).(['stim_' num2str(j) '_correct'])';
        shadedErrorBar(x_ticks,mean(this_traces,1),...
            std(this_traces,[],1),{'color',trial_color.(['correct_stim' num2str(j)]),'linewidth',2},0.1);
        
        this_traces = cell_struct(i).(['stim_' num2str(j) '_incorrect'])';
        shadedErrorBar(x_ticks,mean(this_traces,1),...
            std(this_traces,[],1),{'color',trial_color.(['incorrect_stim' num2str(j)]),'linewidth',2},0.1);
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
plot_cell_type = 'port';
plot_cell_idx = cell_idx_struct.(plot_cell_type);
num_plot_cols = 5;
num_plot_rows = ceil(numel(plot_cell_idx)/num_plot_cols);

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
export_fig([fig_save_path filesep plot_cell_type 'STATrace_' strrep(caiman_file,'.mat','.png')])

%% Plot STAs by imagesc

%% factor analysis
fa_opt.bin_size = 1;
fa_opt.gocue_bin = floor(opt.sta_gocue_frame/fa_opt.bin_size);
fa_opt.stim_bin = ceil(opt.sta_stim_frame/fa_opt.bin_size);
fa_opt.frame_rate = 30/fa_opt.bin_size;
fa_opt.Fs = opt.frame_rate;
fa_opt.trial_length = opt.trial_length/fa_opt.bin_size;
fa_opt.trial_color = trial_color;

fa_opt.idx_fields = {opt.trigger_idx_fd};
fa_opt.fd_names ={'stim_1_correct','stim_1_incorrect','stim_2_correct','stim_2_incorrect'};
fa_opt.plot_fds = fa_opt.fd_names;
fa_opt.m = 4;
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

[fa_traj_struct,num_fs] = get_pop_vectors(fa_struct.F,fa_opt.trial_length,fa_trial_idx);
toc
disp('...Done')
plot_pop_vectors(fa_traj_struct,fa_opt.plot_fds,fa_opt.m,fa_opt,...
    'plot_ylabel','Factor level','plot_num_cols',2);

%% get coding direction
% using raw onine traces to mimic online condition
choice_fds = {'stim_1_correct','stim_2_correct'};
cd_opt = opt;
cd_opt.trial_color = trial_color;
coding_direct = cell2mat(arrayfun(@(x)mean(mean(raw_cell_struct(x).(choice_fds{1})(opt.sta_avg_frames,:),2)),cell_idx_struct.(opt.trigger_idx_fd),'UniformOutput', false))...
    - cell2mat(arrayfun(@(x)mean(mean(raw_cell_struct(x).(choice_fds{2})(opt.sta_avg_frames,:),2)),cell_idx_struct.(opt.trigger_idx_fd),'UniformOutput', false));
% project to coding direction
cd_traj_struct = get_cd_projection(raw_cell_struct(cell_idx_struct.(opt.trigger_idx_fd)),coding_direct);
plot_pop_vectors(cd_traj_struct,choice_fds,1,cd_opt,...
        'plot_ylabel','Projection')
%% choose which to use - factor analysis, or coding direction
if opt.IF_USE_FA
    traj_struct = fa_traj_struct;
else
    traj_struct = cd_traj_struct;
end

%% get stim decoder (using easy trials)
stim_opt = fa_opt;
stim_opt.fd_names = {'stim_1_correct','stim_2_correct'}; % stim 1 will be positive
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
values.stim1 = mean(stim_proj_struct.(stim_opt.fd_names{1})(:,opt.sta_avg_frames),2) ;
values.stim2 = mean(stim_proj_struct.(stim_opt.fd_names{2})(:,opt.sta_avg_frames),2) ;
scatter_cmp_conditions(values,[],...
    1,[0 0 0;0 0 0],'connect_scatter',0,'BriefXlabel',0,'ShowMeanInXlabel',1,'add_jitter',1);
axis square

%% get choice decoder (directly from traj, using threshold trials) 
choice_opt = fa_opt;
choice_opt.fd_names = {'stim_1_correct','stim_1_incorrect'}; %correct choice will be positive
% choice_opt.fd_names = {TS_L2_fds{1},TS_L1_fds{1}};
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

test_fd_names = {'stim_1_correct','stim_1_incorrect',...
                 'stim_2_correct','stim_2_incorrect',...
};

stim_proj_struct = struct();
[stim_proj_struct] = get_projections(raw_cell_struct(cell_idx_struct.(opt.trigger_idx_fd)),stim_norm_weights,test_fd_names,'proj_struct',stim_proj_struct,'bias',-stim_norm_thresh,'IS_CELL_STRUCT',1);

plot_pop_vectors(stim_proj_struct,test_fd_names,1,opt,...
        'ylimit',[-10 10],'plot_ylabel','Projection','plot_num_cols',2,'IF_PLOT_RAW_ONLY',1)
suptitle('Stim decoder projections')
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
    'ylimit',[-10 10],'plot_ylabel','Projection','plot_num_cols',2,'IF_PLOT_RAW_ONLY',1)
suptitle('Choice decoder projections')
%% get choice decoder after projecting to the stim axis - doesn't help with anything
% L2 comes first
% keep IF_CROSSVAL = 1 if sample size is not too small(<5), to avoid biased
% by single trials in the unshuffled field
decod_fds = { 'stim_1_incorrect','stim_1_correct'}; % first fd will be positive; pyRTAOI takes positive stim1 value as incorrect trial
choice_after_stim_struct = struct();
choice_after_stim_struct =  get_binary_classifier( choice_after_stim_struct,stim_proj_struct, choice_opt,...
    'IF_CROSSVAL',1,'IF_FRAMEWISE',0,'fd_names',decod_fds);
choicestim_proj_struct = struct();
 [choicestim_proj_struct] = get_projections(stim_proj_struct,choice_after_stim_struct.B(:,2:end)',decod_fds,'proj_struct',choicestim_proj_struct,'bias',choice_after_stim_struct.B(:,1));
    
plot_pop_vectors(choicestim_proj_struct,decod_fds,1,fa_opt,...
    'ylimit',[-10 10],'plot_ylabel','Projection','plot_num_cols',2,'IF_PLOT_RAW_ONLY',0)

%% save cell structure 
save_brief_fds = {'st_correct_stim_1','st_correct_stim_2','st_incorrect_stim_1','st_incorrect_stim_2','sta_amp_correct_stim_1','sta_amp_correct_stim_2','correct_stimAUC','correct_stimAUC_zscore','is_tuned','cnm_idx','pref_tex','jsf'};
save_brief_fds_names = {'st_correct_stim_1','st_correct_stim_2','st_incorrect_stim_1','st_incorrect_stim_2','sta_amp_correct_stim_1','sta_amp_correct_stim_2','tex_auc','tex_auc_zscore','is_tuned','cnm_idx','pref_tex','jsf'};
save_cell_struct_brief = struct();
for c = 1:num_cells
    for f = 1:numel(save_brief_fds)
        try
        save_cell_struct_brief(c).(save_brief_fds_names{f}) =  cell_struct(c).(save_brief_fds{f});
        catch
            disp([save_brief_fds{f} 'not saved'])
        end
    end
end

tex_output_struct = struct();
tex_output_struct.cell_struct = save_cell_struct_brief;
tex_output_struct.cell_idx_struct = cell_idx_struct;
tex_output_struct.opt = opt;
tex_output_struct.input_caiman_file = fullfile(caiman_path,caiman_file);

output_save_name = [save_path filesep  'ProcTex_' caiman_file];
save(output_save_name,'tex_output_struct')
disp(['Output struct saved as:' output_save_name ])


%% =========================== CHECK PLOTS ======================================
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

