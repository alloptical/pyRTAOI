% generate sta from pyRTAOI result file after a sensory stim movie
% note: rois idx doesn't match with what's shown in pyRTAOI image window
% run this to quickly find sensory-opsin-positve cells  

% ZZ 2018

% adapted from OnlineProcVisual
%% add path - change this for rig
clear all
close all
clc

% ZZ PC
addpath(genpath('C:\Users\User\Desktop\pyRTAOI-rig\Matlab'));
cd('C:\Users\User\Desktop\pyRTAOI-rig\Matlab');

% BRUKER1
% matlab_set_paths_zz

%% stim parameter - CHANGE THIS
sens_duration = 90; % frames
photo_duration = 0; % frames

% params for calculating sta 
opt.N = 1; % threshold for significant auc
opt.all_reset_frames = 30;
opt.pre_exp_frames = 0;
opt.sta_pre_frames = 90;
opt.sta_post_frames = 90;
opt.sta_baseline_frames = 30;
opt.window_size = 60;
opt.sta_avg_frames = 30;
opt.withold_frames_adj = opt.sta_pre_frames:opt.sta_pre_frames+opt.sta_avg_frames;
opt.bs_frame_range = opt.sta_pre_frames+[(-opt.sta_baseline_frames+1):0];
opt.offcell_avg_frames = 30;
opt.sta_thresh = 1;
opt.correct_trial_only = false; % only use correct trials to get tunning
sta_struct = struct();

[trial_color] = online_tex_init_color();
%% load CAIMAN data
[file,path] = uigetfile('*.mat','Select caiman data');
disp(['Loaded file :',fullfile(path,file)])
caiman_data = load(fullfile(path,file)) 
FLAG_PYBEHAV_LOADED = false;
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

%% setup save path and name
opt.output_path = path;
opt.exp_name = strrep(file,'.mat','proc');

%% check if trial order matches in behavior and caiman file
if FLAG_PYBEHAV_LOADED
    if isequal(trials.stim_type,caiman_data.trialOrder)
        disp('data matched')
    else
        warning('trial order mismatch!')
    end
end
%% organise data (generate plots for sanity check)
tex_stim_frames = {};
if(~isempty(caiman_data.stim_frames_caiman))
    sens_stim_frames = caiman_data.sensory_stim_frames+caiman_data.num_frames_init;
else
    sens_stim_frames = [];
end

% trial type
trialOrder = caiman_data.trialOrder;
trialTypes = unique(trialOrder);
num_stim_type = length(unique(trialOrder)); % orientations/texture

% only get correct trials
if opt.correct_trial_only &&  FLAG_PYBEHAV_LOADED
    sens_stim_frames = sens_stim_frames(trials.correct==1);
    trialOrder = trialOrder(trials.correct==1);
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
try
    accepted_idx = caiman_data.accepted+1;
catch
    accepted_idx = caiman_data.accepted_idx+1;
end
opsin_positive_idx = accepted_idx(opsin_positive>0);

% color trial by stim type
hsv_lut = colormap(hsv);
hsv_lut = hsv_lut(2:end-3,:);
indices = round(linspace(1,size(hsv_lut,1),num_stim_type));
opt.trial_color = zeros(numel(trialOrder),3);
for t = 1:numel(trialOrder)
    opt.trial_color(t,:) = trial_color.(['stim' num2str(trialOrder(t))]);
end

num_trials = numel(sens_stim_frames);
opt.tint_factor = 0.6;
opt.type_color = tint(hsv_lut(indices(1:num_stim_type),:),opt.tint_factor);

%% make data structure
cnm_struct = struct();
cnm_dims = caiman_data.cnm_dims;
cnm_image = reshape(caiman_data.cnm_b,cnm_dims);
cnm_A = full(caiman_data.cnm_A);
num_frames = caiman_data.t_cnm;


num_comp = size(cnm_A,2);
comp_shape =[cnm_dims(1),cnm_dims(2),num_comp];
cm = com(sparse(double(cnm_A)),cnm_dims(1),cnm_dims(2)); % center of mass


if ~isempty(caiman_data.frames_skipped)
    skip_frames = caiman_data.frames_skipped + caiman_data.num_frames_init;
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
    
    % putative touch frame (start sta)
    cnm_struct(i).stim_frames = opt.withold_frames_adj(1)+sens_stim_frames;

end

temp_trace = nan(1,tot_frames);
temp_trace(caiman_frames) =  caiman_data.noisyC(1,1:num_frames);
temp_trace = fillmissing(temp_trace,'linear');
backgroundC = temp_trace;

%% only get accepted cells
accepted_idx = caiman_data.accepted_idx+1;
num_cells = numel(accepted_idx);


%% only analyse frames from the cell locked movie
% caiman_data.init_t = 15240+500;
glob_trialon_frames = caiman_data.sensory_stim_frames + caiman_data.t_init;

%% plot spatial components
com_fov = zeros(cnm_dims);
binary_fov = zeros(cnm_dims);
for i = accepted_idx
    com_fov = com_fov+cnm_struct(i).shape;
end

cnm_plot_options = CNMFSetParms;
cnm_plot_options.roi_color = [colormap(lines);colormap(lines);colormap(lines)];

figure('name','fov')
subplot(1,2,1)
imagesc(com_fov)
colormap(gray)
axis square
title('Detected ROIs')


subplot(1,2,2)
[CC,jsf] = plot_contours(sparse(double(cnm_A)),cnm_image,cnm_plot_options,1,[],[],[1 1 1]);
colormap(gray)
axis square
title('GCaMP')

% subplot(1,3,3)
% hold on
% [CC,jsf] = plot_contours(sparse(double(cnm_A)),caiman_data.opsin_mask,cnm_plot_options,1,[],[],[1 1 1]);
% title('Opsin mask')
% set(gca,'YDir','reverse')


% save coords to cell struct
% only take accepted components as cells 
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
    plot(cnm_struct(i).deconvC_full+i*plot_offset,'color','black','linewidth',1.5)
    stim_cell_count = stim_cell_count+1;
end
xlabel('Frames')
ylabel('ROI index')
yticks([ 1:num_comp].*plot_offset)
yticklabels(1:num_comp)
xlim([caiman_data.t_init tot_frames])

% background
plot(backgroundC,'color',[.5 .5 .5],'linestyle',':')

for i = 1:numel(glob_trialon_frames)
    this_color = trial_color.(['stim' num2str(trials.stim_type(i) )]);
    plot([glob_trialon_frames(i) glob_trialon_frames(i)],ylim,'color',this_color)
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
        [~,~,~,~,~,cell_struct(i).sta_traces,cell_struct(i).sta_trace] = make_sta_from_traces(this_cell_trace,this_sens_stim_frames,opt.sta_pre_frames,opt.sta_post_frames,1:opt.sta_baseline_frames);
        cell_struct(i).sta_amp = mean(cell_struct(i).sta_trace(opt.sta_pre_frames:opt.sta_pre_frames+opt.sta_avg_frames));
        
        if  cell_struct(i).sta_amp > opt.sta_thresh
            cell_struct(i).is_sensory = 1;
        end
        
        % get sta for each stim type
        tuning = nan(1,num_stim_type);
        for s = 1:num_stim_type
            [~,~,~,~,~,sta_struct.exp(i).(['type' num2str(s) '_traces']),sta_struct.exp(i).(['type' num2str(s) '_avg'])] =...
                make_sta_from_traces(this_cell_trace,tex_stim_frames{s} ,opt.sta_pre_frames,opt.sta_post_frames,1:opt.sta_baseline_frames);

            % trial average
            tuning(s) = nanmean(sta_struct.exp(i).(['type' num2str(s) '_avg'])(opt.sta_pre_frames:opt.sta_pre_frames+opt.sta_avg_frames));
            cell_struct(i).(['stim', num2str(s) 'amp']) = tuning(s);
            
        end
        

        if max(tuning)>1
            cell_struct(i).is_sensory = 1;
            cell_struct(i).pref_orient = find(tuning == max(tuning));           
        end
        
        cell_struct(i).('tuning') = tuning;
    else
        for s = 1:num_stim_type
             cell_struct(i).(['stim', num2str(s) 'amp']) = nan;
             cell_struct(i).('tuning') = nan(1,num_stim_type);
        end
        
    end
    
end

% sta of backgound component (for alignment check)
[~,~,~,~,~,bg_sta_traces,bg_sta_trace] =...
    make_sta_from_traces(backgroundC,sens_stim_frames ,opt.sta_pre_frames,opt.sta_post_frames,1:opt.sta_baseline_frames);

%% ROC analysis
trial_indices = struct(); % % get trialtype-outcome indices
for i = 1:num_stim_type
    trial_indices.(['correct_stim_' num2str(i)]) = find(trials.correct==1&trials.stim_type==i);
    trial_indices.(['incorrect_stim_' num2str(i)]) = find(trials.incorrect==1&trials.stim_type==i);
    for c = 1:num_cells
        if ~isempty(trial_indices.(['correct_stim_' num2str(i)]))
            cell_struct(c).(['st_correct_stim_' num2str(i)]) = cell_struct(c).sta_traces( trial_indices.(['correct_stim_' num2str(i)]),:)';
        else
            cell_struct(c).(['st_correct_stim_' num2str(i)])  = [];
        end
        if ~isempty(trial_indices.(['incorrect_stim_' num2str(i)]))
            cell_struct(c).(['st_incorrect_stim_' num2str(i)]) = cell_struct(c).sta_traces( trial_indices.(['incorrect_stim_' num2str(i)]),:)';
        else
             cell_struct(c).(['st_incorrect_stim_' num2str(i)]) = [];
        end
    end
end

%% GET CELL IDENTITY
% stimulus AUC calculated from correct trials
num_shuf = 100;
avg_frame_range = 1:30; % temp 
for i = 1:num_cells
all_stim1 = [mean(cell_struct(i).('st_correct_stim_1')(avg_frame_range,:),1)];
all_stim2 = [mean(cell_struct(i).('st_correct_stim_2')(avg_frame_range,:),1)];
labels = [ones(1,length(all_stim1)),2.*ones(1,length(all_stim2))]';
scores = [ all_stim1  all_stim2]';
[~,~,~, correct_stimulusAUC] = perfcurve(labels,scores,1);

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
cell_idx_struct.tex2 = find(correct_stimulusAUC_zscore>opt.N& correct_texAUC>0.6); % cells prefering texture1 in correct trials
cell_idx_struct.tex1 = find(correct_stimulusAUC_zscore<-opt.N& correct_texAUC<0.4); % cells prefering texture2 in correct trials
cell_idx_struct.tex = unique([cell_idx_struct.tex1,cell_idx_struct.tex2]);

for c = 1:num_cells
    cell_struct(i).is_tuned = 0;
    cell_struct(i).pref_tex = 0;
    if ismember(c,cell_idx_struct.tex)
        cell_struct(i).is_tuned = 1;
        if ismember(c,cell_idx_struct.tex1)
            cell_struct(i).pref_tex = 1;
        else
            cell_struct(i).pref_tex = 2;
        end
    end
end

%%
figure('name','correct stim auc')
subplot(1,2,1)
histogram(extractfield(cell_struct,'correct_stimAUC'))

subplot(1,2,2)
histogram(extractfield(cell_struct,'correct_stimAUC_zscore'))



%% MAKE OUTPUT FILE FOR PYRTAOI
opt.target_idx_fd = 'tex';
opt.trigger_idx_fd = 'tex';
opt.fov_size = double(cnm_dims);
opt.ds_factor = caiman_data.ds_factor
[output] = generate_cell_idx_file(cell_struct,cell_idx_struct,[],opt);


%% =========================== PLOTS ======================================
%% Show sensory cells on maps
figure('name','pref. texture on fov');
subplot(1,2,1)
imagesc(com_fov)
colormap(gray)
colorbar('location','southoutside');
axis square
title('Detected ROIs')

ax1 = subplot(1,2,2)
value_field = 'pref_tex';
plot_value_in_rois( cell_struct, value_field,[256 256],ax1,'colorlut',[[1,1,1];opt.type_color],'IF_NORM_PIX',0,'IF_CONTOUR',1,'IF_SHOW_OPSIN',1,'zlimit',[0 4]);
set(gca,'Ydir','reverse')
title('Sensory cells (colored by pref. texture)')

%% Plot STAs for trigger cells
figure('name','trigger sta traces')
num_plot_cols = 4;
num_plot_rois = length(cell_idx_struct.(opt.trigger_idx_fd));
num_plot_rows = ceil(num_plot_rois/num_plot_cols);
plot_count = 1;
for ii = 1:num_plot_rois
    i = cell_idx_struct.(opt.trigger_idx_fd)(ii);
    subtightplot(num_plot_rows,num_plot_cols,plot_count)
    % plot traces
    hold on
    for t = 1:size(cell_struct(i).sta_traces,1)
        plot(cell_struct(i).sta_traces(t,:),'color',opt.trial_color(t,:) ,'linewidth',1)
    end
    plot(cell_struct(i).sta_trace,'color',[.5 .5 .5],'linewidth',1.5)
    
    xlim([0 opt.sta_pre_frames+opt.sta_post_frames+1])
    set(gca,'xtick',[],'xcolor',[1 1 1])
    axis square
    
    plot([opt.sta_pre_frames opt.sta_pre_frames ],ylim,'color',[.5 .5 .5])

    plot_count = plot_count+1;
    
end
suptitle('Trigger cells, stim-triggered response')
%% Show STA on maps
figure('name','sta on fov','position',[200 200 1200 800])
plot_row = 2;
plot_col = ceil(num_stim_type/plot_row);
sta_img = struct();
for s = 1:num_stim_type
    ax1 = subplot(plot_row,plot_col,s);
    value_field = ['stim', num2str(s) 'amp'];
    sta_img.(value_field) = plot_value_in_rois( cell_struct, value_field,[256 256],ax1,'IF_NORM_PIX',0,'IF_CONTOUR',0,'IF_SHOW_OPSIN',1);
    set(gca,'Ydir','reverse')
    title(['Stim' num2str(s)])
end
% export_fig  C:\Users\Zihui\Dropbox\pyRTAOI_demo_figure\visual_stim_inihib\20181029_t016_4sta_on_maps.pdf -painters 

%%
%%%%%%%%%%%%%%%%%%% MORE DETIALED PLOTS GOES BELOW %%%%%%%%%%%%%%%%%%%%%%%%
%% Plot STAs for all components
figure('name','condition sta traces')
num_plot_cols = 6;
num_plot_rows = ceil(num_comp/num_plot_cols);
for i = 1:num_comp
    subtightplot(num_plot_rows,num_plot_cols,i)
    % plot traces
    hold on
    for t = 1:size(cell_struct(i).sta_traces,1)
        plot(cell_struct(i).sta_traces(t,:),'color',opt.trial_color(t,:) ,'linewidth',1)
    end
    plot(cell_struct(i).sta_trace,'color',[.5 .5 .5],'linewidth',1.5)
    
    xlim([0 opt.sta_pre_frames+opt.sta_post_frames+1])
    set(gca,'xtick',[],'xcolor',[1 1 1])
    axis square
    
    if(isempty(find(accepted_idx==i)))
        text(1,1,['ROI ' num2str(i)],'units','normalized','color',[.7 .7 .7],'Horizontalalignment','right','VerticalAlignment','top')
    else
        text(1,1,['ROI ' num2str(i)],'units','normalized','color','black','Horizontalalignment','right','VerticalAlignment','top')
    end
    
    if(~isempty(find(opsin_positive_idx==i)))
        text(1,1,['ROI ' num2str(i)],'units','normalized','color','r','Horizontalalignment','right','VerticalAlignment','top')
    end
    box off
    if( cell_struct(i).is_sensory)
        box on
        set(gca,'XColor',tint([1,0,0],0.5),'YColor',tint([1,0,0],0.5),'linewidth',3)
    end
end
%% Compare targeted and non-targeted sensory cells
sen_target_idx = find([extractfield(cell_struct,'opsin_positive').* extractfield(cell_struct,'accepted').* extractfield(cell_struct,'is_sensory')]>0);
sen_nontarget_idx = find([not(extractfield(cell_struct,'opsin_positive')).* extractfield(cell_struct,'accepted').* extractfield(cell_struct,'is_sensory')]>0);

figure
subplot(1,2,1)
for s = 1:num_stim_type
    this_sta = [];
    for c = sen_target_idx
        if cell_struct(c).pref_orient == s
            this_sta = [this_sta; cell2mat({sta_struct.exp(c).(['type' num2str(cell_struct(c).pref_orient) '_traces'])}')];
        end
    end
    plot(this_sta', 'color',opt.type_color(s,:))
    hold on;
end
plot([opt.sta_baseline_frames opt.sta_baseline_frames],ylim,'color','r')
plot([opt.sta_baseline_frames opt.sta_baseline_frames]+photo_duration,ylim,'color','r')
plot([opt.sta_baseline_frames opt.sta_baseline_frames]+sens_duration,ylim,'color','black')
xlim([0 opt.sta_pre_frames+opt.sta_post_frames+1])
axis square
title('conditioning movie, targets')

subplot(1,2,2)
for s = 1:num_stim_type
    this_sta = [];
    for c = sen_nontarget_idx
        if cell_struct(c).pref_orient == s
            this_sta = [this_sta; cell2mat({sta_struct.exp(c).(['type' num2str(cell_struct(c).pref_orient) '_traces'])}')];
        end
    end
    plot(this_sta', 'color',opt.type_color(s,:))
    hold on;
end

plot([opt.sta_baseline_frames opt.sta_baseline_frames],ylim,'color','r')
plot([opt.sta_baseline_frames opt.sta_baseline_frames]+photo_duration,ylim,'color','r')
plot([opt.sta_baseline_frames opt.sta_baseline_frames]+sens_duration,ylim,'color','black')
xlim([0 opt.sta_pre_frames+opt.sta_post_frames+1])
axis square
title('conditioning movie, non-targets')

%% tuning
target_idx = find([extractfield(cell_struct,'opsin_positive').* extractfield(cell_struct,'accepted')]>0);
nontarget_idx =find( [not(extractfield(cell_struct,'opsin_positive')).* extractfield(cell_struct,'accepted')]>0);

figure
for s = 1:num_stim_type
    subplot(2,ceil(num_stim_type/2),s)
    this_sta = cell2mat({sta_struct.exp(target_idx).(['type' num2str(s) '_avg'])});
    plot(this_sta)
    hold on;
    plot([opt.sta_baseline_frames opt.sta_baseline_frames],ylim,'color','r')
    plot([opt.sta_baseline_frames opt.sta_baseline_frames]+photo_duration,ylim,'color','r')
    plot([opt.sta_baseline_frames opt.sta_baseline_frames]+sens_duration,ylim,'color','black')
    xlim([0 opt.sta_pre_frames+opt.sta_post_frames+1])
    axis square
    title(['Orient.' num2str(s)])
end
suptitle('sensory movie')

%%
figure
hold on
plot(bg_sta_traces')
plot(bg_sta_trace,'color','black','linewidth',2)
 xlim([0 opt.sta_pre_frames+opt.sta_post_frames+1])
 title('back ground sta')
%% plot photo responsive cells on fov
figure
subplot(1,4,1)
imagesc(cnm_image)
colormap(gray)
axis off
colorbar('location','southoutside');
title('GCaMP')
axis square


subplot(1,4,2)
colormap(gray)
imagesc(com_fov)
set(gca,'YDir','normal')
axis off
% plot_contours(sparse(double(cnm_A)),com_fov,cnm_plot_options,1,[],[],[1 1 1]);
colorbar('location','southoutside');
title('Detected ROIs')
axis square


ax1 = subplot(1,4,3);
value_field = 'num_trials';
colorlut = [[1,1,1];opt.trial_color];
plot_value_in_rois( cell_struct, value_field,[256 256],ax1,'colorlut',colorlut,'IF_NORM_PIX',0)
title('Number of visual+photo trials received')

ax2 = subplot(1,4,4);
value_field = 'sta_amp';
colorlut = [];
plot_value_in_rois( cell_struct, value_field,[256 256],ax2,'colorlut',colorlut,'IF_NORM_PIX',0)
title('Visual+Photo')




