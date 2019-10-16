% generate sta from pyRTAOI result file after a sensory stim movie
% note: rois idx doesn't match with what's shown in pyRTAOI image window
% run this to quickly find sensory-opsin-positve cells  

% ZZ 2018

% adapted from OnlineProcVisual

% TO DO:
% concatinate multiple sessions
%
%% add path - change this for rig
clear all
close all
clc

% ZZ PC
% addpath(genpath('C:\Users\User\Desktop\pyRTAOI-rig\Matlab'));
% cd('C:\Users\User\Desktop\pyRTAOI-rig\Matlab');

% BRUKER1
% matlab_set_paths_zz

%% stim parameter - CHANGE THIS
sens_duration = 15; % frames
photo_duration = 0; % frames

opt.N = 1; % threshold for significant auc
opt.sta_pre_frames = 150; 
opt.sta_post_frames = 120;
opt.sta_baseline_frames = 30; % relative to beginning of sta traces

opt.gocue_frame = 135; % relative to trial start
opt.end_rw_frame = 195; % end of response window

% frame indices relative to sta trace
opt.sta_gocue_frame = opt.sta_pre_frames;
opt.sta_trialon_frame = opt.sta_pre_frames-opt.gocue_frame;
opt.sta_avg_frames = [-30:1:0]+opt.sta_gocue_frame;
opt.sta_peak_search_range =  [-30:1:0]+opt.sta_gocue_frame;


opt.withold_frames_adj = [120 150]; 

opt.offcell_avg_frames = 30;
opt.sta_amp_thresh = 1;
opt.frame_rate = 30;

opt.flag_use_peak = true; % if false then will use average to get roc

opt.correct_trial_only = true; % only use correct trials to get tunning

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
sens_stim_frames(sens_stim_frames>caiman_data.t_cnm-opt.sta_post_frames-opt.sta_pre_frames) = [];
num_trials = length(sens_stim_frames);

% trial type
trialOrder = caiman_data.trialOrder(1:num_trials);
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
    cnm_struct(i).noisyC = caiman_data.noisyC(i+1,1:num_frames);
    cnm_struct(i).deconvC = caiman_data.cnm_C(i,1:num_frames);
    cnm_struct(i).centroid = cm(i,:);
    cnm_struct(i).frame_added = find(cnm_struct(i).noisyC >0,1);
    
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

            cell_struct(i).sta_amp = mean(cell_struct(i).sta_trace(opt.sta_avg_frames));
        
        if  cell_struct(i).sta_amp > opt.sta_amp_thresh
            cell_struct(i).is_sensory = 1;
        end        
    end
    
end

% sta of backgound component (for alignment check)
[~,~,~,~,~,bg_sta_traces,bg_sta_trace] =...
    make_sta_from_traces(backgroundC,sens_stim_frames+opt.sta_gocue_frame ,opt.sta_pre_frames,opt.sta_post_frames,1:opt.sta_baseline_frames);

%% ROC analysis
% sta for different trial types
trial_indices = struct(); % % get trialtype-outcome indices
for i = 1:num_stim_type
    trial_indices.(['correct_stim_' num2str(i)]) = find(trials.correct==1&trials.stim_type==i);
    trial_indices.(['incorrect_stim_' num2str(i)]) = find(trials.incorrect==1&trials.stim_type==i);
    for c = 1:num_cells
        if ~isempty(trial_indices.(['correct_stim_' num2str(i)]))
            cell_struct(c).(['st_correct_stim_' num2str(i)]) = cell_struct(c).sta_traces( trial_indices.(['correct_stim_' num2str(i)]),:)';
            cell_struct(c).(['sta_amp_correct_stim_' num2str(i)]) = mean (mean(cell_struct(c).(['st_correct_stim_' num2str(i)])(opt.sta_avg_frames,:),1));

        else
            cell_struct(c).(['st_correct_stim_' num2str(i)])  = [];
            cell_struct(c).(['sta_amp_correct_stim_' num2str(i)]) = [];
        end
        if ~isempty(trial_indices.(['incorrect_stim_' num2str(i)]))
            cell_struct(c).(['st_incorrect_stim_' num2str(i)]) = cell_struct(c).sta_traces( trial_indices.(['incorrect_stim_' num2str(i)]),:)';
            cell_struct(c).(['sta_amp_incorrect_stim_' num2str(i)]) = mean (mean(cell_struct(c).(['st_incorrect_stim_' num2str(i)])(opt.sta_avg_frames,:),1));

        else
             cell_struct(c).(['st_incorrect_stim_' num2str(i)]) = [];
             cell_struct(c).(['sta_amp_incorrect_stim_' num2str(i)]) = 0;
        end
        
    end
end

%% GET CELL IDENTITY
% stimulus AUC calculated from correct trials
num_shuf = 300;
peak_frame_range = opt.sta_peak_search_range;
avg_frame_range = opt.sta_avg_frames;
for i = 1:num_cells
    if ~ opt.flag_use_peak
        all_stim1 = mean(cell_struct(i).('st_correct_stim_1')(avg_frame_range,:),1);
        all_stim2 = mean(cell_struct(i).('st_correct_stim_2')(avg_frame_range,:),1);
    else
        all_stim1 = max(cell_struct(i).('st_correct_stim_1')(peak_frame_range,:),[],1);
        all_stim2 = max(cell_struct(i).('st_correct_stim_2')(peak_frame_range,:),[],1);

    end
labels = [ones(1,length(all_stim1)),2.*ones(1,length(all_stim2))]';
scores = [ all_stim1  all_stim2]';
scores(scores<0.01)=0;
[~,~,~, correct_stimulusAUC] = perfcurve(labels,scores,2);

% shuffle to get zscore stim auc
shuf_stim_auc = nan(1,num_shuf);
for s = 1:num_shuf
    shuf_labels = labels(randperm(length(labels))');
    [~,~,~, shuf_stim_auc(s)] = perfcurve(shuf_labels,scores,1);
end
cell_struct(i).correct_stimAUC_zscore = (correct_stimulusAUC-mean(shuf_stim_auc))/std(shuf_stim_auc);
cell_struct(i).correct_stimAUC = correct_stimulusAUC;
end

correct_stimulusAUC_zscore = extractfield(cell_struct,'correct_stimAUC_zscore');
correct_texAUC = extractfield(cell_struct,'correct_stimAUC');
cell_idx_struct = struct();
cell_idx_struct.tex2 = find(correct_stimulusAUC_zscore>opt.N& correct_texAUC>0.55); % cells prefering texture1 in correct trials
cell_idx_struct.tex1 = find(correct_stimulusAUC_zscore<-opt.N& correct_texAUC<0.45); % cells prefering texture2 in correct trials
cell_idx_struct.tex = unique([cell_idx_struct.tex1,cell_idx_struct.tex2]);

for c = 1:num_cells
    cell_struct(c).is_tuned = 0;
    cell_struct(c).pref_tex = 0;
    if ismember(c,cell_idx_struct.tex)
        cell_struct(c).is_tuned = 1;
        if ismember(c,cell_idx_struct.tex1)
            cell_struct(c).pref_tex = 1;
        else
            cell_struct(c).pref_tex = 2;
        end
    end
end
%% stimulus AUC calculated from incorrect trials
for i = 1:num_cells
    if ~ opt.flag_use_peak
        all_stim1 = mean(cell_struct(i).('st_incorrect_stim_1')(avg_frame_range,:),1);
        all_stim2 = mean(cell_struct(i).('st_incorrect_stim_2')(avg_frame_range,:),1);
    else
        all_stim1 = max(cell_struct(i).('st_incorrect_stim_1')(peak_frame_range,:),[],1);
        all_stim2 = max(cell_struct(i).('st_incorrect_stim_2')(peak_frame_range,:),[],1);

    end
labels = [ones(1,length(all_stim1)),2.*ones(1,length(all_stim2))]';
scores = [ all_stim1  all_stim2]';
scores(scores<0.01)=0;
[~,~,~, incorrect_stimulusAUC] = perfcurve(labels,scores,2);

% shuffle to get zscore stim auc
shuf_stim_auc = nan(1,num_shuf);
for s = 1:num_shuf
    shuf_labels = labels(randperm(length(labels))');
    [~,~,~, shuf_stim_auc(s)] = perfcurve(shuf_labels,scores,1);
end
cell_struct(i).incorrect_stimAUC_zscore = (incorrect_stimulusAUC-mean(shuf_stim_auc))/std(shuf_stim_auc);
cell_struct(i).incorrect_stimAUC = incorrect_stimulusAUC;
end
incorrect_stimulusAUC_zscore = extractfield(cell_struct,'incorrect_stimAUC_zscore');
incorrect_texAUC = extractfield(cell_struct,'incorrect_stimAUC');
cell_idx_struct.port2 = find(correct_stimulusAUC_zscore>opt.N & incorrect_stimulusAUC_zscore<-opt.N); % cells prefering texture1 in correct trials
cell_idx_struct.port1 = find(correct_stimulusAUC_zscore<-opt.N& incorrect_stimulusAUC_zscore>opt.N); % cells prefering texture2 in correct trials
cell_idx_struct.port = unique([cell_idx_struct.port1,cell_idx_struct.port2]);

%% select cell identity for readout and stimulation
opt.target_idx_fd = 'tex';
opt.trigger_idx_fd = 'port';
opt.fov_size = double(cnm_dims);
opt.ds_factor = caiman_data.ds_factor;

%% Plot STAs for all components
figure('name','condition sta traces','units','normalized','outerposition',[0 0 1 1])

num_plot_cols = 6;
num_plot_rows = ceil(num_cells/num_plot_cols);
for i = 1:num_cells
    subtightplot(num_plot_rows,num_plot_cols,i)
    hold on
    % correct trials
    for j = 1:num_stim_type
        plot(cell_struct(i).(['st_correct_stim_' num2str(j)]),'color',trial_color.(['correct_stim' num2str(j)]),'linewidth',1)
        plot(cell_struct(i).(['st_incorrect_stim_' num2str(j)]),'color',trial_color.(['incorrect_stim' num2str(j)]),'linewidth',1)
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
        this_color = trial_color.(['correct_stim' num2str(cell_struct(i).pref_tex)]);
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
export_fig([fig_save_path filesep 'TexSTATrace_' strrep(caiman_file,'.mat','.png')])

%% MAKE OUTPUT FILE FOR PYRTAOI
[output] = generate_cell_idx_file(cell_struct,cell_idx_struct,[],opt);

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

