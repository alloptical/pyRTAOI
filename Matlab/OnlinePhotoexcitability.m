%% process pyrtaoi results after running sequential photostim
% to do: match sensory/photo responsiveness of two groups

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

opt.N = 1; % auc zscore threshold
opt.sta_amp_thresh = 0.1; % sta amplitude threshold

opt.pre_exp_frames = 0;
opt.sta_pre_frames = 20;
opt.sta_post_frames = 30;
opt.sta_baseline_frames = 20;
opt.bs_frame_range = opt.sta_pre_frames+[(-opt.sta_baseline_frames+1):0];
opt.sta_avg_frames = 20; % frames after stim to average for response amp
opt.sta_thresh = 1;
opt.frame_rate = 30; % Hz

cnm_plot_options = CNMFSetParms;
cnm_plot_options.roi_color = [colormap(lines);colormap(lines);colormap(lines)];
close
%% color lut
hsv = colormap(hsv);
hsv = hsv(2:end-3,:);
close 
%% load CAIMAN data
[caiman_file,caiman_path] = uigetfile('*.mat','Select seq photostim caiman data');
disp(['Loaded file :',fullfile(caiman_path,caiman_file)])
caiman_data = load(fullfile(caiman_path,caiman_file));

%% config save path
save_path = [caiman_path filesep 'analysis_files'];
if ~exist(save_path, 'dir')
mkdir(save_path)
end

fig_save_path = [caiman_path filesep 'figures'];
if ~exist(fig_save_path, 'dir')
mkdir(fig_save_path)
end
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

%% load _OutputParams_ mat file 
% in case non target cells are not removed from recording
try
    [file,path] = uigetfile('*.mat','Select cell identity data',caiman_path);
    disp(['Loaded file :',fullfile(path,file)])
    cell_identity_data = load(fullfile(path,file));
    target_cell_idx = cell_identity_data.output.target_idx; % this is idx in cell_struct of OnlineProcTexture
catch
    target_cell_idx = 1:num_cells;
    warning('taking all accepted components as targets')
end

%% make sta based on photostim sequence
opt.frames_with_photo = double(caiman_data.photoDuration/opt.frame_rate+1); % discard skipped frames with photostim when calculating photo sta amplitude
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

[cell_struct] = make_cell_photostim_sta_struct(cnm_struct,cell_struct,accepted_idx,photo_stim_frames,photo_sequence_cell_idx,opt);

for t = 1:numel(target_cell_idx)
    target_cnm_idx(t) = cell_struct(target_cell_idx(t)).cnm_idx; % this is idx in cnm_struct
end

%% make a stim order file for STA movie maker
oris = photo_sequence_cell_idx;
num_diff_stims = length(unique(oris));
save([save_path filesep  'STAMovieMaker_stimType_' num2str(num_diff_stims) '_Stims_' caiman_file],'oris')
% 
%% plot traces
% color trial by stim orientation
unique_stim_type = unique(photo_sequence_cell_idx);
num_stim_type = length(unique_stim_type);
num_stim_trials = length(photo_sequence_cell_idx);
opt.trial_color = zeros(num_stim_trials,3);
opt.tint_factor = 0.6;

indices = round(linspace(1,size(hsv,1),num_stim_type));
opt.type_color = tint(hsv(indices(1:num_stim_type),:),opt.tint_factor);
for s = 1:length(unique_stim_type)
    this_stim_type = unique_stim_type(s);
    opt.trial_color(photo_sequence_cell_idx == this_stim_type,:) = repmat(opt.type_color(s,:) ,numel(find(photo_sequence_cell_idx==this_stim_type)),1);
end

figure; hold on
plot_offset = 5;
stim_cell_count = 1;
non_stim_cell_count = 1;

for i = 1:num_comp
    plot(zscore(cnm_struct(i).deconvC_full)+i*plot_offset,'color',[.7 .7 .7],'linewidth',1.5)
    stim_cell_count = stim_cell_count+1;
end

for i = target_cnm_idx
    plot(zscore(cnm_struct(i).deconvC_full)+i*plot_offset,'color','black','linewidth',1.5)
    stim_cell_count = stim_cell_count+1;
end

xlabel('Frames')
ylabel('ROI index')
yticks([ 1:num_comp].*plot_offset)
yticklabels(1:num_comp)

for i = 1:num_stim_trials
    plot([photo_stim_frames(i) photo_stim_frames(i)],ylim,'color',opt.trial_color(i,:),'linewidth',1)
end

for s = 1:num_stim_type
    plot([caiman_data.t_init tot_frames],target_cnm_idx(s)*plot_offset.*[1 1],'color',opt.type_color(s,:),'linewidth',.5)
end
xlim([caiman_data.t_init tot_frames])
ylim([0 num_comp].*plot_offset+3)
%% plot fov with sta amp
figure('name','fov','units','normalized','outerposition',[0 0 1 1])
[CC,jsf] = plot_contours(sparse(double(full(caiman_data.cnm_A))),cnm_image,cnm_plot_options,1,[],[],[1 1 1]);
colormap(gray)
axis square
title('GCaMP')
close
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
    cell_struct(i).shape = cnm_struct(this_idx).shape;
    cell_struct(i).opsin_positive = 0;
    cell_struct(i).cnm_idx = this_idx;
    cell_struct(i).jsf = jsf(this_idx);
end
%% plot spatial components and photostim sta
target_com_fov = zeros(cnm_dims);
binary_fov = zeros(cnm_dims);
for i = 1: numel(target_cell_idx)
    target_com_fov = target_com_fov+cell_struct(target_cell_idx(i)).shape;
end

cnm_plot_options = CNMFSetParms;
cnm_plot_options.roi_color = [colormap(lines);colormap(lines);colormap(lines)];

figure('name','sta fov','units','normalized','outerposition',[0 0 1 1])
subplot(1,4,1)
imagesc(com_fov)
colormap(gray)
colorbar('location','southoutside');
axis square
title('Detected ROIs')

subplot(1,4,2)
plot_contours(sparse(double(full(caiman_data.cnm_A))),cnm_image,cnm_plot_options,1,[],[],[1 1 1]);
colormap(gray)
colorbar('location','southoutside');
axis square
title('GCaMP')

subplot(1,4,3)
imagesc(target_com_fov)
colormap(gray)
colorbar('location','southoutside');
axis square
title('Targeted ROIs')


ax1 = subplot(1,4,4);
plot_value_in_rois( cell_struct, 'sta_amp',[256 256],ax1,'IF_NORM_PIX',0,'IF_CONTOUR',0,'IF_SHOW_OPSIN',0,'show_cell_idx',target_cell_idx);
set(gca,'Ydir','reverse')
title('Target photo response')

export_fig([fig_save_path filesep 'PhotoSTAFOV_' strrep(caiman_file,'.mat','.png')])

%% plot sta traces
figure('name','target sta traces','units','normalized','outerposition',[0 0 1 1])
num_plot_cols = 5;
num_plot_rois = length(target_cell_idx);
num_plot_rows = ceil(num_plot_rois/num_plot_cols);
plot_count = 1;
for ii = 1:num_plot_rois
    i = target_cell_idx(plot_count);
    subtightplot(num_plot_rows,num_plot_cols,plot_count)
    % plot traces
    hold on
    for t = 1:size(cell_struct(i).sta_traces,1)
        plot(cell_struct(i).sta_traces(t,:),'color',[.5 .5 .5] ,'linewidth',1)
    end
    plot(cell_struct(i).sta_trace,'color',[.5 .5 .5],'linewidth',1.5)
    
    xlim([0 opt.sta_pre_frames+opt.sta_post_frames+1])
    set(gca,'xtick',[],'xcolor',[1 1 1])
    axis square
    
    plot([opt.sta_pre_frames opt.sta_pre_frames ],ylim,'color',[1,0,0],'linestyle',':','linewidth',2) % photostim start frame
    plot([1,1].* opt.sta_pre_frames+opt.frames_with_photo,ylim,'linestyle',':','color',[1,0,0]) % end of photostim
    plot(opt.frames_with_photo+opt.sta_pre_frames + [1,opt.sta_avg_frames],[0 0],'color',[0,0,0],'linewidth',2) % window used to calculate photostim amplitude

    
    
    plot_count = plot_count+1;  
    if cell_struct(i).is_photo
        text_color = [1,0,0];
    else
        text_color = [0 0 0];
    end
    text(0.05,1,['Cell ' num2str(i) ' (ROI ' num2str(cell_struct(i).cnm_idx) ')'],'units','normalized', 'horizontalalignment','left', 'color',text_color)
    text(0.05,.9,['photostim auc ' num2str(cell_struct(i).photo_auc,'%10.2f')],'units','normalized', 'horizontalalignment','left', 'color',text_color)
    text(0.05,.8,['zscore auc '  num2str(cell_struct(i).photo_auc_zscore,'%10.2f') ],'units','normalized', 'horizontalalignment','left', 'color',text_color)
        
end

export_fig([fig_save_path filesep 'PhotoSTATrace' strrep(caiman_file,'.mat','.png')])

%% save structure
% make a brief struct to combine with texture struct later
save_brief_fds = {'sta_amp','photo_auc','photo_auc_zscore','is_photo','cnm_idx','jsf'};
save_brief_fds_names = {'photo_sta_amp','photo_auc','photo_auc_zscore','is_photo','cnm_idx','jsf'};
save_cell_struct_brief = struct();
for c = 1:num_cells
    for f = 1:numel(save_brief_fds)
        save_cell_struct_brief(c).(save_brief_fds_names{f}) =  cell_struct(c).(save_brief_fds{f});
    end
end
photo_output_struct = struct();
photo_output_struct.cell_struct = save_cell_struct_brief;
photo_output_struct.opt = opt;
photo_output_struct.input_caiman_file = fullfile(caiman_path,caiman_file);

output_save_name = [save_path filesep  'ProcPhotoExcit_' caiman_file];
save(output_save_name,'photo_output_struct')
disp(['Output struct saved as:' output_save_name ])