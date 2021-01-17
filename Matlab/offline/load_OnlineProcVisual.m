% load data
caiman_data = load('Y:\zzhang\Data\20181024\pyrtaoi_results\20181024_L527_t_0011_rtaoi_DS_2.0_OnlineProc_DS_2.0_165202.mat') % - example

%% color lut
hsv = colormap(hsv);
hsv = hsv(2:end-3,:);
close 

%% stim types
num_stim_type = 8; % orientations
num_stim_per_type = 4;
vis_duration = 30; % frames
photo_duration = 90; % frames
vis_stim_frames = {};
init_vis_stim_frames = {};

if(~isempty(caiman_data.stim_frames_caiman))
    sens_stim_frames = caiman_data.sensory_stim_frames+caiman_data.num_frames_init;
else
    sens_stim_frames = [];
end
for i = 1:num_stim_type
    vis_stim_frames{i} = sens_stim_frames(i-1+[1:num_stim_type:num_stim_type*num_stim_per_type]);
    init_vis_stim_frames{i} = vis_stim_frames{i} - caiman_data.num_frames_init;
end
%% indices
opsin_positive = caiman_data.opsin_positive;
try
    accepted_idx = caiman_data.accepted+1;
catch
    accepted_idx = caiman_data.accepted_idx+1;
end
opsin_positive_idx = accepted_idx(opsin_positive>0);
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



%% plot spatial components
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
hold on
[CC,jsf] = plot_contours(sparse(double(cnm_A)),caiman_data.opsin_mask,cnm_plot_options,1,[],[],[1 1 1]);
title('Opsin mask')
set(gca,'YDir','reverse')

subplot(1,3,3)
plot_contours(sparse(double(cnm_A)),cnm_image,cnm_plot_options,1,[],[],[1 1 1]);



%% save coords to cell struct
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


%% plot traces
figure; hold on
plot_offset = 5;
stim_cell_count = 1;
non_stim_cell_count = 1;

opt.trial_color = [];
num_trials = numel(sens_stim_frames);
opt.trial_color = zeros(num_trials,3);
opt.tint_factor = 0.6;

indices = round(linspace(1,size(hsv,1),num_stim_type));

% color trial by stim orientation
opt.type_color = tint(hsv(indices(1:num_stim_type),:),opt.tint_factor);
for s = 1:num_stim_type
opt.trial_color(s+[1:num_stim_type:num_stim_type*num_stim_per_type]-1,:) = repmat(opt.type_color(s,:) ,num_stim_per_type,1);
end


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
plot((num_comp+1)*plot_offset + backgroundC,'color',[.5 .5 .5])


for i = 1:numel(glob_stim_frames)
    plot([glob_stim_frames(i) glob_stim_frames(i)],ylim,'color',opt.trial_color(i,:))
end

for i = 1:numel(init_stim_frames)
    plot([init_stim_frames(i) init_stim_frames(i)],ylim,'color',opt.trial_color(i,:))
end


%% stim triggered average 
opt.N = 2;
opt.all_reset_frames = 30;
opt.pre_exp_frames = 0;
opt.opt.window_size = 60;
opt.sta_pre_frames = 30;
opt.sta_post_frames = 120;
opt.sta_baseline_frames = 30;
opt.window_size = 60;
opt.sta_avg_frames = 30;
opt.offcell_avg_frames = 30;
opt.sta_thresh = 1;
sta_struct = struct();

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
    
    if(~isempty(find(accepted_idx==i)))
        cell_struct(i).accepted = 1;
    end
    
    if(this_num_trials>0)
        % average across all stim types
        [~,~,~,~,~,cell_struct(i).sta_traces,cell_struct(i).sta_trace] = make_sta_from_traces(this_cell_trace,sens_stim_frames,opt.sta_pre_frames,opt.sta_post_frames,1:opt.sta_baseline_frames);
        cell_struct(i).sta_amp = mean(cell_struct(i).sta_trace(opt.sta_pre_frames:opt.sta_pre_frames+opt.sta_avg_frames));
        cell_struct(i).off_amp = mean(cell_struct(i).sta_trace(opt.sta_pre_frames+vis_duration+[1:opt.offcell_avg_frames]));
        
        [~,~,~,~,~,cell_struct(i).bs_sta_traces,cell_struct(i).bs_sta_trace] = make_sta_from_traces(this_cell_trace,init_stim_frames,opt.sta_pre_frames,opt.sta_post_frames,1:opt.sta_baseline_frames);
        cell_struct(i).bs_sta_amp = mean(cell_struct(i).bs_sta_trace(opt.sta_pre_frames:opt.sta_pre_frames+opt.sta_avg_frames));
        cell_struct(i).bs_off_amp = mean(cell_struct(i).bs_sta_trace(opt.sta_pre_frames+vis_duration+[1:opt.offcell_avg_frames]));
        
        if  cell_struct(i).sta_amp > opt.sta_thresh
            cell_struct(i).is_sensory = 1;
        end
              
        % get sta for each stim type
        tuning_ctr = nan(1,num_stim_type);
        tuning_exp = tuning_ctr;
        for s = 1:num_stim_type
            [~,~,~,~,~,sta_struct.exp(i).(['type' num2str(s) '_traces']),sta_struct.exp(i).(['type' num2str(s) '_avg'])] =...
                make_sta_from_traces(this_cell_trace,vis_stim_frames{s} ,opt.sta_pre_frames,opt.sta_post_frames,1:opt.sta_baseline_frames);
            tuning_exp(s) = nanmean(sta_struct.exp(i).(['type' num2str(s) '_avg'])(opt.sta_pre_frames:opt.sta_pre_frames+opt.sta_avg_frames));

            [~,~,~,~,~,sta_struct.ctr(i).(['type' num2str(s) '_traces']),sta_struct.ctr(i).(['type' num2str(s) '_avg'])] = ...
                make_sta_from_traces(this_cell_trace,init_vis_stim_frames{s} ,opt.sta_pre_frames,opt.sta_post_frames,1:opt.sta_baseline_frames);
            tuning_ctr(s) = nanmean(sta_struct.ctr(i).(['type' num2str(s) '_avg'])(opt.sta_pre_frames:opt.sta_pre_frames+opt.sta_avg_frames));
            cell_struct(i).(['stim', num2str(s) 'amp']) = tuning_exp(s);
            cell_struct(i).(['bs_stim', num2str(s) 'amp']) = tuning_ctr(s);
           
        end
        
        if max(tuning_ctr)>1
            cell_struct(i).is_sensory = 1;
            cell_struct(i).pref_orient = find(tuning_ctr == max(tuning_ctr));
        end
        
        if (cell_struct(i).bs_off_amp>1) && (~ cell_struct(i).is_sensory) 
            cell_struct(i).is_offcell = 1;
        end
        cell_struct(i).('bs_tuning') = tuning_ctr;
        cell_struct(i).('tuning') = tuning_exp;
    end
    
end

% sta of backgound component (for alignment check)
[~,~,~,~,~,bg_sta_traces,bg_sta_trace] =...
    make_sta_from_traces(backgroundC,sens_stim_frames ,opt.sta_pre_frames,opt.sta_post_frames,1:opt.sta_baseline_frames);

%% Show STA on maps
figure('name','baseline sta on fov','position',[200 200 1200 800])
plot_row = 2;
plot_col = ceil(num_stim_type/plot_row);
for s = 1:num_stim_type
    ax1 = subplot(plot_row,plot_col,s)
    value_field = ['bs_stim', num2str(s) 'amp'];
    plot_value_in_rois( cell_struct, value_field,[256 256],ax1,'IF_NORM_PIX',0,'IF_CONTOUR',1,'IF_SHOW_OPSIN',1)
    title(['Stim' num2str(s)])
end

figure('name','conditioning sta on fov','position',[200 200 1200 800])
plot_row = 2;
plot_col = ceil(num_stim_type/plot_row);
for s = 1:num_stim_type
    ax1 = subplot(plot_row,plot_col,s)
    value_field = ['stim', num2str(s) 'amp'];
    plot_value_in_rois( cell_struct, value_field,[256 256],ax1,'IF_NORM_PIX',0,'IF_CONTOUR',1,'IF_SHOW_OPSIN',1)
    title(['Stim' num2str(s)])
end

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
    
    xlim([0 length(cell_struct(i).sta_trace)])
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
figure('name','baseline sta traces')
num_plot_cols = 6;
num_plot_rows = ceil(num_comp/num_plot_cols);
for i = 1:num_comp
    subtightplot(num_plot_rows,num_plot_cols,i)
    % plot traces
    hold on
    for t = 1:size(cell_struct(i).bs_sta_traces,1)
        plot(cell_struct(i).bs_sta_traces(t,:),'color',opt.trial_color(t,:) ,'linewidth',1)
    end
    plot(cell_struct(i).bs_sta_trace,'color',[.5 .5 .5],'linewidth',1.5)
    
    xlim([0 length(cell_struct(i).sta_trace)])
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
subplot(2,2,1)
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
xlim([0 length(cell_struct(i).sta_trace)])
axis square
title('conditioning movie, targets')

subplot(2,2,2)
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
xlim([0 length(cell_struct(i).sta_trace)])
axis square
title('conditioning movie, non-targets')

subplot(2,2,3)
for s = 1:num_stim_type
    this_sta = [];
    for c = sen_target_idx
        if cell_struct(c).pref_orient == s
            this_sta = [this_sta; cell2mat({sta_struct.ctr(c).(['type' num2str(cell_struct(c).pref_orient) '_traces'])}')];
        end
    end
    plot(this_sta', 'color',opt.type_color(s,:))
    hold on;
end
plot([opt.sta_baseline_frames opt.sta_baseline_frames],ylim,'color','r')
plot([opt.sta_baseline_frames opt.sta_baseline_frames]+photo_duration,ylim,'color','r')
plot([opt.sta_baseline_frames opt.sta_baseline_frames]+vis_duration,ylim,'color','black')
xlim([0 length(cell_struct(i).sta_trace)])
axis square
title('initialisation movie, targets')

subplot(2,2,4)
for s = 1:num_stim_type
    this_sta = [];
    for c = sen_nontarget_idx
        if cell_struct(c).pref_orient == s
            this_sta = [this_sta; cell2mat({sta_struct.ctr(c).(['type' num2str(cell_struct(c).pref_orient) '_traces'])}')];
        end
    end
    plot(this_sta', 'color',opt.type_color(s,:))
    hold on;
end
plot([opt.sta_baseline_frames opt.sta_baseline_frames],ylim,'color','r')
plot([opt.sta_baseline_frames opt.sta_baseline_frames]+photo_duration,ylim,'color','r')
plot([opt.sta_baseline_frames opt.sta_baseline_frames]+vis_duration,ylim,'color','black')
xlim([0 length(cell_struct(i).sta_trace)])
axis square
title('initialisation movie, non-targets')


%%
figure
for s = 1:num_stim_type
    subplot(2,ceil(num_stim_type/2),s)
    this_sta = cell2mat({sta_struct.exp(sen_target_idx).(['type' num2str(s) '_avg'])});
    plot(this_sta)
    hold on;
    plot([opt.sta_baseline_frames opt.sta_baseline_frames],ylim,'color','r')
    plot([opt.sta_baseline_frames opt.sta_baseline_frames]+photo_duration,ylim,'color','r')
    plot([opt.sta_baseline_frames opt.sta_baseline_frames]+vis_duration,ylim,'color','black')
    xlim([0 length(cell_struct(i).sta_trace)])
    axis square
    title(['Orient.' num2str(s)])

end
suptitle('conditioning movie, targets')

figure
for s = 1:num_stim_type
    subplot(2,ceil(num_stim_type/2),s)
    this_sta = cell2mat({sta_struct.exp(sen_nontarget_idx).(['type' num2str(s) '_avg'])});
    plot(this_sta)
    hold on;
    plot([opt.sta_baseline_frames opt.sta_baseline_frames],ylim,'color','r')
    plot([opt.sta_baseline_frames opt.sta_baseline_frames]+photo_duration,ylim,'color','r')
    plot([opt.sta_baseline_frames opt.sta_baseline_frames]+vis_duration,ylim,'color','black')
    xlim([0 length(cell_struct(i).sta_trace)])
    axis square
    title(['Orient.' num2str(s)])
end
suptitle('conditioning movie, non-targets')

figure
for s = 1:num_stim_type
    subplot(2,ceil(num_stim_type/2),s)
    this_sta = cell2mat({sta_struct.ctr(sen_target_idx).(['type' num2str(s) '_avg'])});
    plot(this_sta)
    hold on;
    plot([opt.sta_baseline_frames opt.sta_baseline_frames],ylim,'color','r')
    plot([opt.sta_baseline_frames opt.sta_baseline_frames]+photo_duration,ylim,'color','r')
    plot([opt.sta_baseline_frames opt.sta_baseline_frames]+vis_duration,ylim,'color','black')
    xlim([0 length(cell_struct(i).sta_trace)])
    axis square
    title(['Orient.' num2str(s)])
end
suptitle('initialisation movie, targets')

figure
for s = 1:num_stim_type
    subplot(2,ceil(num_stim_type/2),s)
    this_sta = cell2mat({sta_struct.ctr(sen_nontarget_idx).(['type' num2str(s) '_avg'])});
    plot(this_sta)
    hold on;
    plot([opt.sta_baseline_frames opt.sta_baseline_frames],ylim,'color','r')
    plot([opt.sta_baseline_frames opt.sta_baseline_frames]+photo_duration,ylim,'color','r')
    plot([opt.sta_baseline_frames opt.sta_baseline_frames]+vis_duration,ylim,'color','black')
    xlim([0 length(cell_struct(i).sta_trace)])
    axis square
    title(['Orient.' num2str(s)])
end
suptitle('initialisation movie, non-targets')

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
    xlim([0 length(cell_struct(i).sta_trace)])
    axis square
    title(['Orient.' num2str(s)])
end
suptitle('conditioning movie')

figure
for s = 1:num_stim_type
    subplot(2,ceil(num_stim_type/2),s)
    this_sta = cell2mat({sta_struct.ctr(target_idx).(['type' num2str(s) '_avg'])});
    plot(this_sta)
    hold on;
    plot([opt.sta_baseline_frames opt.sta_baseline_frames],ylim,'color','r')
    plot([opt.sta_baseline_frames opt.sta_baseline_frames]+photo_duration,ylim,'color','r')
    plot([opt.sta_baseline_frames opt.sta_baseline_frames]+vis_duration,ylim,'color','black')
    xlim([0 length(cell_struct(i).sta_trace)])
    axis square
    title(['Orient.' num2str(s)])
end
suptitle('initialisation movie')
%%
figure
hold on
plot(bg_sta_traces')
plot(bg_sta_trace,'color','black','linewidth',2)
xlim([0 length(cell_struct(i).sta_trace)])
title('back ground sta')
%%
figure;
hold on
subplot(1,2,1)
title('baseline tuning')
ctr_mat = cell2mat({cell_struct.bs_tuning}');
exp_mat = cell2mat({cell_struct.tuning}');
climit = [min([ctr_mat(:);exp_mat(:)]),max([ctr_mat(:);exp_mat(:)])];
imagesc(ctr_mat,climit)
colorbar('location','southoutside');
subplot(1,2,2)
title('+photostim tuning')
imagesc(exp_mat,climit)
colorbar('location','southoutside');


%% plot photo responsive cells on fov
figure
subplot(1,5,1)
imagesc(cnm_image)
colormap(gray)
axis off
colorbar('location','southoutside');
title('GCaMP')
axis square


subplot(1,5,2)
colormap(gray)
imagesc(com_fov)
set(gca,'YDir','normal')
axis off
% plot_contours(sparse(double(cnm_A)),com_fov,cnm_plot_options,1,[],[],[1 1 1]);
colorbar('location','southoutside');
title('Detected ROIs')
axis square


ax1 = subplot(1,5,3);
value_field = 'num_trials';
colorlut = [[1,1,1];opt.trial_color];
plot_value_in_rois( cell_struct, value_field,[256 256],ax1,'colorlut',colorlut,'IF_NORM_PIX',0)
title('Number of visual+photo trials received')

ax2 = subplot(1,5,4);
value_field = 'sta_amp';
colorlut = [];
plot_value_in_rois( cell_struct, value_field,[256 256],ax2,'colorlut',colorlut,'IF_NORM_PIX',0)
title('Visual+Photo')

ax3 = subplot(1,5,5);
value_field = 'bs_sta_amp';
colorlut = [];
plot_value_in_rois( cell_struct, value_field,[256 256],ax3,'colorlut',colorlut,'IF_NORM_PIX',0)
title('Visual alone')



