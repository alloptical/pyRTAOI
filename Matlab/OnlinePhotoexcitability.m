%% process pyrtaoi results after running sequential photostim
% to do: match sensory/photo responsiveness of two groups

%% add path - change this for rig
clear all
close all
clc

% ZZ PC
addpath(genpath('C:\Users\User\Desktop\pyRTAOI-rig\Matlab'));
cd('C:\Users\User\Desktop\pyRTAOI-rig\Matlab');

% BRUKER1
% matlab_set_paths_zz

%% parameters - CHANGE THIS
% params for calculating sta 
opt.N = 1;
opt.pre_exp_frames = 0;
opt.sta_pre_frames = 30;
opt.sta_post_frames = 60;
opt.sta_baseline_frames = 30;
opt.bs_frame_range = opt.sta_pre_frames+[(-opt.sta_baseline_frames+1):0];
opt.sta_avg_frames = 15; % frames after stim to average for response amp
opt.sta_thresh = 1;

cnm_plot_options = CNMFSetParms;
cnm_plot_options.roi_color = [colormap(lines);colormap(lines);colormap(lines)];

%% load CAIMAN data
[file,path] = uigetfile('*.mat','Select caiman data');
disp(['Loaded file :',fullfile(path,file)])
caiman_data = load(fullfile(path,file)) 

%% make data struct
[cnm_struct,cnm_image,num_comp,cnm_dims,tot_frames] = make_cnm_struct(caiman_data);
% spatial components
com_fov = zeros(cnm_dims);
binary_fov = zeros(cnm_dims);
for i = 1:num_comp
    com_fov = com_fov+cnm_struct(i).shape;
end

%% load target cell idx file
% in case non target cells are not removed from recording
try
    [file,path] = uigetfile('*.mat','Select cell identity data');
    disp(['Loaded file :',fullfile(path,file)])
    cell_identity_data = load(fullfile(path,file));
    target_idx = cell_identity_data.output.target_idx;
catch
    target_idx = 1:num_comp;
end

%% make sta based on photostim sequence
accepted_idx = caiman_data.accepted_idx+1;
num_cells = numel(accepted_idx);

cell_struct = struct();
photo_stim_frames =  caiman_data.online_photo_frames + caiman_data.t_init;
temp_photo_sequence_idx = caiman_data.photo_sequence_idx+1;
photo_sequence_idx = nan(1,length(photo_stim_frames));
for t = 1:numel(target_idx)
    photo_sequence_idx(temp_photo_sequence_idx == t) = target_idx(t);
end

[cell_struct] = make_cell_photostim_sta_struct(cnm_struct,cell_struct,accepted_idx,photo_stim_frames,photo_sequence_idx,opt);

%% find photoresponsive cells

%% plot traces
figure; hold on
plot_offset = 5;
stim_cell_count = 1;
non_stim_cell_count = 1;

for i = 1:num_comp
    plot(cnm_struct(i).deconvC_full+i*plot_offset,'color',[.5 .5 .5],'linewidth',1.5)
    stim_cell_count = stim_cell_count+1;
end

for i = target_idx
    plot(cnm_struct(i).deconvC_full+i*plot_offset,'color','black','linewidth',1.5)
    stim_cell_count = stim_cell_count+1;
end

xlabel('Frames')
ylabel('ROI index')
yticks([ 1:num_comp].*plot_offset)
yticklabels(1:num_comp)
xlim([caiman_data.t_init tot_frames])

for i = 1:numel(photo_stim_frames)
    plot([photo_stim_frames(i) photo_stim_frames(i)],ylim,'color',[1 0 0],'linestyle',':')
end


%% plot fov with sta amp
figure('name','fov')
subplot(1,2,1)
imagesc(com_fov)
colormap(gray)
axis square
title('Detected ROIs')

subplot(1,2,2)
[CC,jsf] = plot_contours(sparse(double(full(caiman_data.cnm_A))),cnm_image,cnm_plot_options,1,[],[],[1 1 1]);
colormap(gray)
axis square
title('GCaMP')

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
    if(~isempty(find(opsin_positive_idx==i)))
         cell_struct(i).opsin_positive = 1;
    end
end


%% plot spatial components
target_com_fov = zeros(cnm_dims);
binary_fov = zeros(cnm_dims);
for i = 1: numel(target_idx)
    target_com_fov = target_com_fov+cell_struct(target_idx(i)).shape;
end

cnm_plot_options = CNMFSetParms;
cnm_plot_options.roi_color = [colormap(lines);colormap(lines);colormap(lines)];

figure('name','fov')
subplot(1,3,1)
imagesc(com_fov)
colormap(gray)
colorbar('location','southoutside');
axis square
title('Detected ROIs')


subplot(1,3,2)
imagesc(target_com_fov)
colormap(gray)
colorbar('location','southoutside');
axis square
title('Targeted ROIs')

subplot(1,3,3)
ax1 = subplot(1,3,3);
plot_value_in_rois( cell_struct, 'sta_amp',[256 256],ax1,'IF_NORM_PIX',0,'IF_CONTOUR',0,'IF_SHOW_OPSIN',0,'jsf',jsf);
set(gca,'Ydir','reverse')
title('Target photo response')


%%
