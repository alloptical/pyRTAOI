% load data
% caiman_data = load('D:\pyRTAOI data\stim_at_fixed_frames\20180811_OG300_t_0004_DS_2.0_rtaoi_OnlineProc_163822.mat')
% caiman_data = load('D:\pyRTAOI data\stim_at_fixed_frames\20180822_OG299_t_0002_rtaoi_DS_2.0.tirtaoi_OnlineProc_DS_2.0_170940.mat')
% caiman_data = load('D:\pyRTAOI data\stim_at_fixed_frames\20180830_OG323_t_0004_rtaoi_DS_2.0_OnlineProc_DS_2.0_163232.mat')
% caiman_data = load('D:\pyRTAOI data\stim_at_fixed_frames\20180907_OG328_t_0007_rtaoi_DS_2.0_OnlineProc_DS_2.0_160113.mat')
caiman_data = load('D:\pyRTAOI_data\stim_at_fixed_frames\GCaMP6f\20180909_OG347_t_0003_rtaoi_DS_2.0_OnlineProc_DS_2.0_140615.mat') % - example_
% caiman_data = load('D:\pyRTAOI data\stim_at_fixed_frames\GCaMP6f\20180910_OG349_t_0004_rtaoi_DS_2.0_OnlineProc_DS_2.0_163941.mat') 
% caiman_data = load('D:\pyRTAOI data\stim_at_fixed_frames\20180910_OG349_t_0005_rtaoi_DS_2.0_OnlineProc_DS_2.0_164300.mat') 

%% color lut
hsv = colormap(hsv);
hsv = hsv(2:end-3,:);
% close all
%% make data structure
cnm_struct = struct();
cnm_dims = caiman_data.cnm_dims;
cnm_image = reshape(caiman_data.cnm_b,cnm_dims);
cnm_A = full(caiman_data.cnm_A);
num_frames = caiman_data.t_cnm;

num_comp = size(cnm_A,2);
comp_shape =[cnm_dims(1),cnm_dims(2),num_comp];
cm = com(sparse(double(cnm_A)),cnm_dims(1),cnm_dims(2)); % center of mass

if(~isempty(caiman_data.photo_stim_frames_caiman))
    stim_frames = caiman_data.photo_stim_frames_caiman+caiman_data.num_frames_init;
else
    stim_frames = [];
end


for i = 1:num_comp
    cnm_struct(i).shape = reshape(cnm_A(:,i),cnm_dims);
    cnm_struct(i).noisyC = caiman_data.noisyC(i+1,1:num_frames);
    cnm_struct(i).deconvC = caiman_data.cnm_C(i,1:num_frames);
    cnm_struct(i).centroid = cm(i,:);
    cnm_struct(i).frame_added = find(cnm_struct(i).noisyC >0,1);
    if(~isempty(stim_frames))
        cnm_struct(i).stim_frames = stim_frames(stim_frames>cnm_struct(i).frame_added);
    else
        cnm_struct(i).stim_frames = [];
    end
end

%% indices
opsin_positive = caiman_data.opsin_positive;
try
    accepted_idx = caiman_data.accepted+1;
catch
    accepted_idx = caiman_data.accepted_idx+1;
end
opsin_positive_idx = accepted_idx(opsin_positive>0);

%% plot spatial components
com_fov = zeros(cnm_dims);
binary_fov = zeros(cnm_dims);
for i = 1:num_comp
    com_fov = com_fov+cnm_struct(i).shape;
end

cnm_plot_options = CNMFSetParms;
cnm_plot_options.roi_color = [colormap(lines);colormap(lines)];

figure('name','fov')
subplot(1,3,1)
imagesc(com_fov)
colormap(gray)
axis square
title('Detected ROIs')

subplot(1,3,2)
hold on
[CC,jsf] = plot_contours(sparse(double(cnm_A)),caiman_data.opsin_mask,cnm_plot_options,1,[],[],[1 1 1]);
title('Opsin mask')
set(gca,'YDir','reverse')

subplot(1,3,3)
plot_contours(sparse(double(cnm_A)),cnm_image,cnm_plot_options,1,[],[],[1 1 1]);

figure
colormap(gray)
plot_contours(sparse(double(cnm_A)),com_fov,cnm_plot_options,1,[],[],[1 1 1]);
axis square

%% save coords to cell struct
cell_struct = struct();
for i = 1:size(jsf,1)
    cell_struct(i).coordinates = jsf(i).coordinates;
    cell_struct(i).pix_values = jsf(i).values;
    cell_struct(i).centroid = jsf(i).centroid;
    
    cell_struct(i).opsin_positive = 0;
    if(~isempty(find(opsin_positive_idx==i)))
         cell_struct(i).opsin_positive = 1;
    end
end


%% plot traces
figure; hold on
plot_offset = 5;
stim_cell_count = 1;
non_stim_cell_count = 1;
sub1_i = 0;
sub2_i = 0;
for i = 1:num_comp
    % cells not stimulated
    if cell_struct(i).opsin_positive == 0
        subplot(3,1,1)
        hold on
%         plot(cnm_struct(i).noisyC+i*plot_offset,'color',cnm_plot_options.roi_color(i,:))
        plot(cnm_struct(i).deconvC+sub1_i*plot_offset,'color',[.5 .5 .5],'linewidth',1.5)
        non_stim_cell_count = non_stim_cell_count+1;
        sub1_i = sub1_i+1;
    else
        subplot(3,1,2:3)
        hold on
%         plot(cnm_struct(i).noisyC+i*plot_offset,'color',cnm_plot_options.roi_color(i,:))
        plot(cnm_struct(i).deconvC+sub2_i*plot_offset,'color','black','linewidth',1.5)
        stim_cell_count = stim_cell_count+1;
        sub2_i = sub2_i+1;
    end
end
xlabel('Frames')
ylabel('ROI index')
yticks([ 1:num_comp].*plot_offset)
yticklabels(1:num_comp)
xlim([0 size(cnm_struct(1).noisyC,2)])

for i = 1:numel(stim_frames)
    plot([stim_frames(i) stim_frames(i)],ylim,'r')
end
% export_fig  C:\Users\Zihui\Dropbox\pyRTAOI_demo_figure\photostim_fix_frames\20180909_t003_traces.pdf -painters 

%% plot traces by imagesc

opsin_positive = find(extractfield(cell_struct,'opsin_positive')==1);
opsin_cell_traces = cell2mat({cnm_struct(opsin_positive).deconvC}');
ds_factor = 10;
ds_opsin_cell_traces = nan(size(opsin_cell_traces,1),ceil(size(opsin_cell_traces,2)/ds_factor));
% downsample for print
for r = 1:size(opsin_cell_traces,1)
ds_opsin_cell_traces(r,:) = downsample(opsin_cell_traces(r,:),10);
end

figure('position',[100 100 1200 300])
imagesc(ds_opsin_cell_traces)
colormap(gray)
colormap(b2r(min(ds_opsin_cell_traces(:)),max(ds_opsin_cell_traces(:))))
colorbar('eastoutside')
set(gca,'YDir','reverse')
hold on
axis off
for i = 1:numel(stim_frames)
    plot([stim_frames(i) stim_frames(i)],ylim,'black')
end
% export_fig  C:\Users\Zihui\Dropbox\pyRTAOI_demo_figure\photostim_fix_frames\20180909_t003_traces_red.pdf -painters 


%% stim triggered average
opt.N = 2;
opt.all_reset_frames = 30;
opt.pre_exp_frames = 0;
opt.opt.window_size = 60;
opt.sta_pre_frames = 30;
opt.sta_post_frames = 60;
opt.sta_baseline_frames = 30;
opt.window_size = 60;
opt.sta_avg_frames = 30;
opt.sta_thresh = 1;
opt.trial_color = [];
num_trials = numel(stim_frames);
opt.trial_color = zeros(num_trials,3);
opt.tint_factor = 0.7;

indices = round(linspace(1,size(hsv,1),num_trials));
for i = 1:num_trials
opt.trial_color(i,:) = tint(hsv(indices(i),:),opt.tint_factor);
end


figure('name','sta traces')
num_plot_cols = 6;
num_plot_rows = ceil(num_comp/num_plot_cols);
for i = 1:num_comp
    subtightplot(num_plot_rows,num_plot_cols,i)
    this_cell_trace = cnm_struct(i).deconvC;
    this_num_trials = numel(cnm_struct(i).stim_frames );
    cell_struct(i).num_trials = this_num_trials*cell_struct(i).opsin_positive;
    cell_struct(i).is_photo = 0;
    cell_struct(i).sta_amp = 0;
    cell_struct(i).sta_traces = [];
    cell_struct(i).sta_trace = [];
    
    if(this_num_trials>0)
        [~,~,~,~,~,cell_struct(i).sta_traces,cell_struct(i).sta_trace] = make_sta_from_traces(this_cell_trace,cnm_struct(i).stim_frames ,opt.sta_pre_frames,opt.sta_post_frames,1:opt.sta_baseline_frames);
        cell_struct(i).sta_amp = mean(cell_struct(i).sta_trace(opt.sta_pre_frames:opt.sta_pre_frames+opt.sta_avg_frames));
        if  cell_struct(i).sta_amp > opt.sta_thresh
            cell_struct(i).is_photo = 1;
        end
        
        % plot traces
        hold on
        for t = 1:this_num_trials
            plot(cell_struct(i).sta_traces(t,:),'color',opt.trial_color(t,:) ,'linewidth',1)
        end
        plot(cell_struct(i).sta_trace,'color',[.5 .5 .5],'linewidth',1.5)
    end
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
    if( cell_struct(i).is_photo )
        box on
        set(gca,'XColor',tint([1,0,0],0.5),'YColor',tint([1,0,0],0.5),'linewidth',3)        
    end
    
end


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
title('Number of photostims received')

ax2 = subplot(1,4,4);
value_field = 'sta_amp';
colorlut = [];
plot_value_in_rois( cell_struct, value_field,[256 256],ax2,'colorlut',colorlut,'IF_NORM_PIX',0)
title('Photostim-triggered average')

suptitle('20180909 OG347 GCaMP6f+C1V1')

% export_fig  D:\pyRTAOI_data\stim_at_fixed_frames\GCaMP6f\plots\20180909_t003_FOV.pdf -painters 


%%
ax2 = figure
value_field = 'sta_amp';
colorlut = [];
plot_value_in_rois( cell_struct, value_field,[256 256],ax2,'colorlut',colorlut,'IF_NORM_PIX',0, 'zlimit',[-3 20])
title('Photostim-triggered average')


