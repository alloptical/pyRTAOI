function [ cell_struct ] = proc_OnlineProcFixedFrames( caiman_data, save_name,opt)
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
    cnm_struct(i).onlineC = caiman_data.online_C(i+1,1:num_frames);
    cnm_struct(i).deconvC = caiman_data.cnm_C(i,1:num_frames);
    cnm_struct(i).centroid = cm(i,:);
    cnm_struct(i).frame_added = find(cnm_struct(i).onlineC >0,1);
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
for i = 1:num_comp
    com_fov = com_fov+cnm_struct(i).shape;
end

cnm_plot_options = CNMFSetParms;
cnm_plot_options.roi_color = colormap(lines);

figure('name','fov')
subplot(1,3,1)
imagesc(com_fov)
colormap(gray)
axis square
title('Detected ROIs')

try
    subplot(1,3,2)
    hold on
    plot_contours(sparse(double(cnm_A)),caiman_data.opsin_mask,cnm_plot_options,1,[],[],[1 1 1]);
    title('Opsin mask')
    set(gca,'YDir','reverse')
catch
    warning('c1v1 mask not found')
end

subplot(1,3,3)
[CC,jsf] = plot_contours(sparse(double(cnm_A)),cnm_image,cnm_plot_options,1,[],[],[1 1 1]);
title('Detected ROIs')

% figure
% colormap(gray)
% plot_contours(sparse(double(cnm_A)),com_fov,cnm_plot_options,1,[],[],[1 1 1]);
% axis square

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
%% stim triggered average
figure('position',[50 50 1200 900],'name','sta traces')
num_plot_cols = 8;
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

print([opt.fig_save_path filesep save_name '_sta.png'],'-dpng');

%% summary plot: fov and traces
figure('name','traces','position',[50 50 1200 900]); hold on
plot_offset = 5;
stim_cell_count = 1;
non_stim_cell_count = 1;

subplot(2,4,1)
imagesc(cnm_image)
colormap(gray)
axis off
colorbar('location','southoutside');
title('GCaMP')
axis square


subplot(2,4,2)
colormap(gray)
imagesc(com_fov)
set(gca,'YDir','normal')
axis off
% plot_contours(sparse(double(cnm_A)),com_fov,cnm_plot_options,1,[],[],[1 1 1]);
colorbar('location','southoutside');
title('Detected ROIs')
axis square


ax1 = subplot(2,4,3);
value_field = 'num_trials';
colorlut = [[1,1,1];opt.trial_color];
plot_value_in_rois( cell_struct, value_field,[256 256],ax1,'colorlut',colorlut,'IF_NORM_PIX',0)
title('Number of photostims received')

ax2 = subplot(2,4,4);
value_field = 'sta_amp';
colorlut = [];
plot_value_in_rois( cell_struct, value_field,[256 256],ax2,'colorlut',colorlut,'IF_NORM_PIX',0)
title('Photostim-triggered average')


subplot(2,4,5:8)
hold on
for i = 1:num_comp
    % cells not stimulated
    if cell_struct(i).opsin_positive == 0
        plot(cnm_struct(i).onlineC+i*plot_offset,'color',cnm_plot_options.roi_color(i,:))
        plot(cnm_struct(i).deconvC+i*plot_offset,'color',[.5 .5 .5],'linewidth',1.5)
        non_stim_cell_count = non_stim_cell_count+1;
    else
        plot(cnm_struct(i).onlineC+i*plot_offset,'color',cnm_plot_options.roi_color(i,:))
        plot(cnm_struct(i).deconvC+i*plot_offset,'color','black','linewidth',1.5)
        stim_cell_count = stim_cell_count+1;

    end
end
xlabel('Frames')
ylabel('ROI index')
yticks([ 1:num_comp].*plot_offset)
yticklabels(1:num_comp)
xlim([0 size(cnm_struct(1).onlineC,2)])

for i = 1:numel(stim_frames)
    plot([stim_frames(i) stim_frames(i)],ylim,'r')
end


suptitle(strrep(save_name,'_',' '))
print([opt.fig_save_path filesep save_name '_summary.png'],'-dpng');

end

