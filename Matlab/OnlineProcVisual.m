% generate sta from pyRTAOI result file after a sensory stim movie
% note: rois idx doesn't match with what's shown in pyRTAOI image window
% run this to quickly find sensory-opsin-positve cells  

% ZZ 2018
%% add path - change this for rig
clear all
close all
clc
addpath(genpath('C:\Users\User\Desktop\pyRTAOI-rig\Matlab'));
cd('C:\Users\User\Desktop\pyRTAOI-rig\Matlab');

%%
% load data
[file,path] = uigetfile('*.mat');
disp(['Loaded file :',fullfile(path,file)])
caiman_data = load(fullfile(path,file)) 


%% stim parameter - CHANGE THIS
num_stim_type = 4; % orientations
num_stim_per_type = 10;
vis_duration = 30; % frames
photo_duration = 90; % frames

% params for calculating sta 
opt.N = 2;
opt.all_reset_frames = 30;
opt.pre_exp_frames = 0;
opt.sta_pre_frames = 30;
opt.sta_post_frames = 120;
opt.sta_baseline_frames = 30;
opt.window_size = 60;
opt.sta_avg_frames = 30;
opt.offcell_avg_frames = 30;
opt.sta_thresh = 1;
sta_struct = struct();

%% organise data (with generate plots for sanity check)
vis_stim_frames = {};
if(~isempty(caiman_data.stim_frames_caiman))
    sens_stim_frames = caiman_data.sensory_stim_frames+caiman_data.num_frames_init;
else
    sens_stim_frames = [];
end
for i = 1:num_stim_type
    vis_stim_frames{i} = sens_stim_frames(i-1+[1:num_stim_type:num_stim_type*num_stim_per_type]);
end

% indices
opsin_positive = caiman_data.opsin_positive;
try
    accepted_idx = caiman_data.accepted+1;
catch
    accepted_idx = caiman_data.accepted_idx+1;
end
opsin_positive_idx = accepted_idx(opsin_positive>0);

% color trial by stim orientation
hsv_lut = colormap(hsv);
hsv_lut = hsv_lut(2:end-3,:);
indices = round(linspace(1,size(hsv_lut,1),num_stim_type));
opt.trial_color = [];
num_trials = numel(sens_stim_frames);
opt.trial_color = zeros(num_trials,3);
opt.tint_factor = 0.6;
opt.type_color = tint(hsv_lut(indices(1:num_stim_type),:),opt.tint_factor);
for s = 1:num_stim_type
opt.trial_color(s+[1:num_stim_type:num_stim_type*num_stim_per_type]-1,:) = repmat(opt.type_color(s,:) ,num_stim_per_type,1);
end

% make data structure
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
glob_stim_frames = caiman_data.sensory_stim_frames + caiman_data.num_frames_init;
init_stim_frames = caiman_data.sensory_stim_frames; % movie used for initialisation had same sensory stim

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
    
    % caiman stim frames
    if(~isempty(sens_stim_frames))
        cnm_struct(i).stim_frames = sens_stim_frames(sens_stim_frames>cnm_struct(i).frame_added);
    else
        cnm_struct(i).stim_frames = [];
    end
end

temp_trace = nan(1,tot_frames);
temp_trace(caiman_frames) =  caiman_data.noisyC(1,1:num_frames);
temp_trace = fillmissing(temp_trace,'linear');
backgroundC = temp_trace;


% plot spatial components
com_fov = zeros(cnm_dims);
binary_fov = zeros(cnm_dims);
for i = 1:num_comp
    com_fov = com_fov+cnm_struct(i).shape;
end

cnm_plot_options = CNMFSetParms;
cnm_plot_options.roi_color = [colormap(lines);colormap(lines);colormap(lines)];

figure('name','fov')
subplot(1,3,1)
imagesc(com_fov)
colormap(gray)
axis square
title('Detected ROIs')


subplot(1,3,2)
plot_contours(sparse(double(cnm_A)),cnm_image,cnm_plot_options,1,[],[],[1 1 1]);
colormap(gray)
axis square
title('GCaMP')

subplot(1,3,3)
hold on
[CC,jsf] = plot_contours(sparse(double(cnm_A)),caiman_data.opsin_mask,cnm_plot_options,1,[],[],[1 1 1]);
title('Opsin mask')
set(gca,'YDir','reverse')


% save coords to cell struct
cell_struct = struct();
for i = 1:size(jsf,1)
    temp_coords = jsf(i).coordinates;
    lin_idx = zeros(size(temp_coords,1),1);
   
    for t = 1:size(temp_coords,1)
        lin_idx(t) = sub2ind(cnm_dims,temp_coords(t,1),temp_coords(t,2));
    end
    
    cell_struct(i).contour = CC{i};
    cell_struct(i).lin_coords = lin_idx;
    cell_struct(i).coordinates = jsf(i).coordinates;
    cell_struct(i).pix_values = jsf(i).values;
    cell_struct(i).centroid = jsf(i).centroid;
    cell_struct(i).opsin_positive = 0;
    if(~isempty(find(opsin_positive_idx==i)))
         cell_struct(i).opsin_positive = 1;
    end
end

% plot traces
figure; hold on
plot_offset = 5;
stim_cell_count = 1;
non_stim_cell_count = 1;

for i = 1:num_comp
    % cells not stimulated
    if cell_struct(i).opsin_positive == 0
%         subplot(2,1,1)
%         hold on
%         plot(cnm_struct(i).noisyC_full+i*plot_offset,'color',cnm_plot_options.roi_color(i,:))
        plot(cnm_struct(i).deconvC_full+i*plot_offset,'color',[.5 .5 .5],'linewidth',1.5)
        non_stim_cell_count = non_stim_cell_count+1;
    else
%         subplot(2,1,2)
%         hold on
%         plot(cnm_struct(i).noisyC_full+i*plot_offset,'color',cnm_plot_options.roi_color(i,:))
        plot(cnm_struct(i).deconvC_full+i*plot_offset,'color','black','linewidth',1.5)
        stim_cell_count = stim_cell_count+1;

    end
end
xlabel('Frames')
ylabel('ROI index')
yticks([ 1:num_comp].*plot_offset)
yticklabels(1:num_comp)
xlim([0 tot_frames])

% background
plot(backgroundC,'color',[.5 .5 .5],'linestyle',':')

for i = 1:numel(glob_stim_frames)
    plot([glob_stim_frames(i) glob_stim_frames(i)],ylim,'color',opt.trial_color(i,:))
end

for i = 1:numel(glob_stim_frames)
    plot([glob_stim_frames(i) glob_stim_frames(i)]+vis_duration,ylim,'color',opt.trial_color(i,:),'linestyle',':')
end
% export_fig  C:\Users\Zihui\Dropbox\pyRTAOI_demo_figure\visual_stim_inihib\20181029_t016_traces_blacklines.pdf -painters 

%% traces as imagesc - comment this out, for poster only
all_traces = cell2mat({cnm_struct.deconvC_full}');
ds_factor = 10;
ds_all_traces = nan(size(all_traces,1),ceil(size(all_traces,2)/ds_factor));
% downsample for print
for r = 1:size(all_traces,1)
ds_all_traces(r,:) = downsample(all_traces(r,:),10);
end

figure('position',[100 100 1200 300])
imagesc(ds_all_traces)
colormap(gray)
colormap(b2r(0,35))
colorbar('eastoutside')
set(gca,'YDir','reverse')
hold on
axis off

for i = 1:numel(glob_stim_frames)
    plot([glob_stim_frames(i) glob_stim_frames(i)]./ds_factor,ylim,'color',opt.trial_color(i,:))
end

for i = 1:numel(glob_stim_frames)
    plot([[glob_stim_frames(i) glob_stim_frames(i)]+vis_duration]./ds_factor,ylim,'color',opt.trial_color(i,:),'linestyle',':')
end
% export_fig  C:\Users\Zihui\Dropbox\pyRTAOI_demo_figure\visual_stim_inihib\20181029_t016_traces_red.pdf -painters 


%% stim triggered average 
for i = 1:num_comp
    this_cell_trace = cnm_struct(i).deconvC_full;
    this_num_trials = numel(cnm_struct(i).stim_frames );
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
        [~,~,~,~,~,cell_struct(i).sta_traces,cell_struct(i).sta_trace] = make_sta_from_traces(this_cell_trace,sens_stim_frames,opt.sta_pre_frames,opt.sta_post_frames,1:opt.sta_baseline_frames);
        cell_struct(i).sta_amp = mean(cell_struct(i).sta_trace(opt.sta_pre_frames:opt.sta_pre_frames+opt.sta_avg_frames));
        
        if  cell_struct(i).sta_amp > opt.sta_thresh
            cell_struct(i).is_sensory = 1;
        end
        
        % get sta for each stim type
        tuning = nan(1,num_stim_type);
        offtuning = nan(1,num_stim_type);
        for s = 1:num_stim_type
            [~,~,~,~,~,sta_struct.exp(i).(['type' num2str(s) '_traces']),sta_struct.exp(i).(['type' num2str(s) '_avg'])] =...
                make_sta_from_traces(this_cell_trace,vis_stim_frames{s} ,opt.sta_pre_frames,opt.sta_post_frames,1:opt.sta_baseline_frames);
            tuning(s) = nanmean(sta_struct.exp(i).(['type' num2str(s) '_avg'])(opt.sta_pre_frames:opt.sta_pre_frames+opt.sta_avg_frames));
            offtuning(s) = nanmean(sta_struct.exp(i).(['type' num2str(s) '_avg'])(opt.sta_pre_frames+vis_duration+[1:opt.offcell_avg_frames]));
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

%% Show sensory cells on maps
figure('name','pref. orientation on fov');

subplot(1,3,1)
imagesc(com_fov)
colormap(gray)
colorbar('location','southoutside');
axis square
title('Detected ROIs')


subplot(1,3,2)
% plot_contours(sparse(double(cnm_A)),caiman_data.opsin_mask,cnm_plot_options,0,[],[],[1 1 1]);
imagesc(caiman_data.opsin_mask)
axis square
colorbar('location','southoutside');
set(gca,'YDir','reverse')
title('Opsin mask')

ax1 = subplot(1,3,3)
value_field = 'pref_orient';
plot_value_in_rois( cell_struct, value_field,[256 256],ax1,'colorlut',[[1,1,1];opt.type_color],'IF_NORM_PIX',0,'IF_CONTOUR',1,'IF_SHOW_OPSIN',1,'zlimit',[0 4]);
set(gca,'Ydir','reverse')
title('Sensory cells (colored by pref. orientation)')

% export_fig  C:\Users\Zihui\Dropbox\pyRTAOI_demo_figure\visual_stim_inihib\20181029_t016_FOV_maps.pdf -painters 
% %% example -delete later
% ax1 = figure
% temp_colors = [[251 185 47];[169 207 137];[147 109 173];[242 151 189]]./255;
% value_field = 'pref_orient';
% plot_value_in_rois( cell_struct, value_field,[256 256],ax1,'colorlut',[[1,1,1];temp_colors],'IF_NORM_PIX',0,'IF_CONTOUR',1,'IF_SHOW_OPSIN',1,'zlimit',[0 4]);
% set(gca,'Ydir','reverse')
% title('Sensory cells (colored by pref. orientation)')
% %  export_fig  C:\Users\Zihui\Dropbox\pyRTAOI_demo_figure\visual_stim_inihib\20181029_t016_FOV_uglycolored.pdf -painters 
%% Plot STAs for visual cells
figure('name','condition sta traces')
num_plot_cols = 4;
num_sensory = length(find(extractfield(cell_struct,'is_sensory')==1));
num_plot_rows = ceil(num_sensory/num_plot_cols);
plot_count = 1;
for i = 1:num_comp
    if cell_struct(i).is_sensory == 0
        continue
    end
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
    plot_count = plot_count+1;
    
end
suptitle('Sensory cells, stim-triggered response')
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
plot([opt.sta_baseline_frames opt.sta_baseline_frames]+vis_duration,ylim,'color','black')
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
plot([opt.sta_baseline_frames opt.sta_baseline_frames]+vis_duration,ylim,'color','black')
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
    plot([opt.sta_baseline_frames opt.sta_baseline_frames]+vis_duration,ylim,'color','black')
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




