% load data

%  caiman_data = load('D:\pyRTAOI data\replay\20180912_OG347_t_0013_rtaoi_DS_2.0_OnlineProc_DS_2.0_182518.mat') 
% caiman_data = load('D:\pyRTAOI data\replay\20180912_OG347_t_0014_rtaoi_DS_2.0_OnlineProc_DS_2.0_183159.mat') 
% caiman_data.num_frames_init = 1348; % 20180912_OG347_t004 - used t013 results; 

% caiman_data = load('D:\pyRTAOI test data\Data\pyrtaoi_results\20180912_OG347_t-014_rtaoi_OnlineProc_DS_2.0_113038.mat');
% caiman_data = load('D:\pyRTAOI data\replay\20180926_L506_t_0002_rtaoi_DS_2.0_OnlineProc_DS_2.0_191152.mat'); % visual


% caiman_data = load('D:\pyRTAOI test data\Data\20180926_L506_t_0010_rtaoi_DS_2.0_OnlineProc_DS_2.0_194609.mat'); % visual
caiman_data = load('D:\pyRTAOI test data\Data\20180926_L506_t_0009_rtaoi_DS_2.0_OnlineProc_DS_2.0_194100.mat'); % visual

sen_stim_duration = 30; % frames
photo_stim_duration = 10;

%% color lut
hsv = colormap(hsv);
hsv = hsv(2:end-3,:);
% close all
%% make data structure
cell_struct = struct();
cnm_dims = caiman_data.cnm_dims;
cnm_image = reshape(caiman_data.cnm_b,cnm_dims);
cnm_A = full(caiman_data.cnm_A);
num_frames = caiman_data.t_cnm;

num_comp = size(cnm_A,2);
comp_shape =[cnm_dims(1),cnm_dims(2),num_comp];
cm = com(sparse(double(cnm_A)),cnm_dims(1),cnm_dims(2)); % center of mass

% photostim
if(~isempty(caiman_data.photo_stim_frames_caiman))
    photo_stim_frames = caiman_data.photo_stim_frames_caiman+caiman_data.num_frames_init;
else
    photo_stim_frames = [];
end

% sensory stim
if(~isempty(caiman_data.stim_frames_caiman))
    sen_stim_frames = caiman_data.stim_frames_caiman+caiman_data.num_frames_init;
else
    sen_stim_frames = [];
end

online_photo_targets = [];
try
    online_photo_targets = cellfun(@(x)(x+1),caiman_data.online_photo_targets,'UniformOutput',false);
catch
    warning('no online photostim targets found')
end

for i = 1:num_comp
    cell_struct(i).shape = reshape(cnm_A(:,i),cnm_dims);
    cell_struct(i).noisyC = caiman_data.noisyC(i+1,1:num_frames); % first row is background
    cell_struct(i).deconvC = caiman_data.cnm_C(i,1:num_frames);
    cell_struct(i).onlineC = caiman_data.online_C(i+1,1:num_frames); % first row is background
    cell_struct(i).centroid = cm(i,:);
    cell_struct(i).frame_added = find(cell_struct(i).noisyC >0,1);
end

%% indices
opsin_positive = caiman_data.opsin_positive;
try
    accepted_idx = caiman_data.accepted+1;
catch
    accepted_idx = caiman_data.accepted_idx+1;
end
opsin_positive_idx = accepted_idx(opsin_positive>0);
num_accepted = numel(accepted_idx);
%% plot spatial components
com_fov = zeros(cnm_dims);
binary_fov = zeros(cnm_dims);
for i = 1:num_comp
    com_fov = com_fov+cell_struct(i).shape;
end

cnm_plot_options = CNMFSetParms;
cnm_plot_options.roi_color = [colormap(lines);colormap(lines)];

figure
colormap(gray)
[CC,jsf] = plot_contours(sparse(double(cnm_A)),com_fov,cnm_plot_options,1,[],[],[1 1 1]);
axis square

figure('name','fov')
subplot(1,2,1)
imagesc(com_fov)
colormap(gray)
axis square

% try
% subplot(1,3,2)
% hold on
% plot_contours(sparse(double(cnm_A)),caiman_data.opsin_mask,cnm_plot_options,1,[],[],[1 1 1]);
% title('Opsin mask')
% set(gca,'YDir','reverse')
% catch
%     warning('opsin mask plot error')
% end

subplot(1,2,2)
plot_contours(sparse(double(cnm_A)),cnm_image,cnm_plot_options,1,[],[],[1 1 1]);

suptitle('Detected ROIs')

%% threshold online
all_thresh = zeros(num_accepted,num_frames);
FLAG_THRESH_FOUND = 0;
try
    for f = 1:num_frames-caiman_data.num_frames_init
        this_values = caiman_data.online_thresh{f};
        this_num_roi = numel(this_values);
        all_thresh(1:this_num_roi,f+caiman_data.num_frames_init) = this_values;
        
    end
    FLAG_THRESH_FOUND = 1;
catch
    warning('no online threshold found')
end

%% save coords to cell struct
inter_trial_frames = max(diff(sen_stim_frames));
for i = 1:num_comp
    cell_struct(i).coordinates = jsf(i).coordinates;
    cell_struct(i).pix_values = jsf(i).values;
    cell_struct(i).centroid = jsf(i).centroid;
    cell_struct(i).online_idx = find(accepted_idx==i);
    
    if(~isempty(cell_struct(i).online_idx))
        cell_struct(i).online_thresh = all_thresh(cell_struct(i).online_idx,:);
        
    end
    
    cell_struct(i).opsin_positive = 0;
    if(~isempty(find(opsin_positive_idx==i)))
         cell_struct(i).opsin_positive = 1;
    end
    
    if(~isempty(online_photo_targets))
        cell_struct(i).photostim_frames = photo_stim_frames(cell2mat(cellfun(@(x)ismember(cell_struct(i).online_idx,x),online_photo_targets,'UniformOutput',false)));
    else
        cell_struct(i).photostim_frames  = [];
    end
    
    cell_struct(i).sen_stim_frames = sen_stim_frames(sen_stim_frames>cell_struct(i).frame_added )
    
    
    temp_idx1 = cell2mat(arrayfun(@(x)find(cell_struct(i).sen_stim_frames<x, 1, 'last' ),cell_struct(i).photostim_frames,'UniformOutput',false));
    temp_idx2 = cell2mat(arrayfun(@(x)find(cell_struct(i).sen_stim_frames>x-inter_trial_frames, 1, 'first' ),cell_struct(i).photostim_frames,'UniformOutput',false));
    
    try
        cell_struct(i).trig_sen_stim_frames = cell_struct(i).sen_stim_frames(intersect(temp_idx1,temp_idx2));
    catch
        cell_struct(i).trig_sen_stim_frames = [];
    end
end


%% plot traces
% sen_stim_frames = [830 1424 2018 2612 3206 3800 4394 4988 5582 6176] % for test

figure('name','ca traces'); hold on
plot_offset = 5;
stim_cell_count = 1;
non_stim_cell_count = 1;
ylimit = [0 (num_comp+1)*plot_offset];

C_field = 'onlineC';

for i = 1:numel(sen_stim_frames)
    plot([sen_stim_frames(i) sen_stim_frames(i)],ylimit,'c','linewidth',2)
    plot([sen_stim_frames(i) sen_stim_frames(i)]+sen_stim_duration,ylimit,'c')
end


for i = 1:num_comp
    % cells not stimulated
    if isempty(cell_struct(i).online_idx)
%         subplot(2,1,1)
%         hold on
%          plot(cell_struct(i).onlineC+i*plot_offset,'color',cnm_plot_options.roi_color(i,:))
         plot(cell_struct(i).(C_field)+i*plot_offset,'color',[.5 .5 .5],'linewidth',1.5)
        non_stim_cell_count = non_stim_cell_count+1;
    else
%         subplot(2,1,2)
%         hold on
%          plot(cell_struct(i).onlineC+i*plot_offset,'color',cnm_plot_options.roi_color(i,:))
         plot(cell_struct(i).(C_field)+i*plot_offset,'color','black','linewidth',1.5)
        stim_cell_count = stim_cell_count+1;

    end
    
    % photostim frames
    if ~isempty(cell_struct(i).photostim_frames)
        for f = 1:numel(cell_struct(i).photostim_frames)
            this_frame = cell_struct(i).photostim_frames(f);
            plot([this_frame this_frame],[i*plot_offset,i*plot_offset+plot_offset/2],'r','linewidth',2)
        end
    end
    
    % thresh
%      plot(cell_struct(i).online_thresh+i*plot_offset,'color','black','linewidth',0.5,'linestyle','--')
     
     % online events
     if((~isempty(cell_struct(i).online_thresh))&&FLAG_THRESH_FOUND)
         event_trace = nan(size(cell_struct(i).(C_field)));
         event_trace(cell_struct(i).(C_field)>cell_struct(i).online_thresh) = cell_struct(i).(C_field)(cell_struct(i).(C_field)>cell_struct(i).online_thresh);
         plot(event_trace+i*plot_offset,'color','r','linewidth',2)
         
         
         % mark trials with events detected
         for s = 1:numel(sen_stim_frames)
             this_window = [sen_stim_frames(s) sen_stim_frames(s)+sen_stim_duration];
             if(ismember(0,isnan(event_trace(this_window(1):this_window(2)))))
                 plot(this_window,[i*plot_offset,i*plot_offset],'m','linewidth',2)
             end
             
         end
     end
     
end
xlabel('Frames')
ylabel('ROI index')
yticks([1:num_comp].*plot_offset)
yticklabels(1:num_comp)
xlim([0 size(cell_struct(1).noisyC,2)])

% for i = 1:numel(photo_stim_frames)
%     plot([photo_stim_frames(i) photo_stim_frames(i)],ylim,'r')
% end

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
num_trials = max([numel(photo_stim_frames) numel(sen_stim_frames)]);
opt.trial_color = zeros(num_trials,3);
opt.tint_factor = 0.7;

indices = round(linspace(1,size(hsv,1),num_trials));
for i = 1:num_trials
opt.trial_color(i,:) = tint(hsv(indices(i),:),opt.tint_factor);
end

figure('name','photo sta traces');

num_plot_cols = 6;
num_plot_rows = ceil(num_comp/num_plot_cols);
for i = 1:num_comp
    subtightplot(num_plot_rows,num_plot_cols,i)
    this_cell_trace = cell_struct(i).deconvC;
    this_num_trials = numel(cell_struct(i).photostim_frames);
    cell_struct(i).num_trials = this_num_trials;
    cell_struct(i).is_photo = 0;
    cell_struct(i).is_sensory = 0;
    cell_struct(i).sta_amp = 0;
    cell_struct(i).sta_traces = [];
    cell_struct(i).sta_trace = [];
    
    % sensory trials that triggered photostim
    if(~isempty(cell_struct(i).trig_sen_stim_frames))
        [~,~,~,~,~,cell_struct(i).trig_sen_sta_traces,cell_struct(i).trig_sen_sta_trace] = make_sta_from_traces(this_cell_trace,cell_struct(i).trig_sen_stim_frames,opt.sta_pre_frames,opt.sta_post_frames,1:opt.sta_baseline_frames);
        cell_struct(i).trig_sens_sta_amp = mean(cell_struct(i).trig_sen_sta_trace(opt.sta_pre_frames:opt.sta_pre_frames+opt.sta_avg_frames));
        
    end
    
    % sensory stim sta 
    [~,~,~,~,~,cell_struct(i).sens_sta_traces,cell_struct(i).sens_sta_trace] = ...
        make_sta_from_traces(this_cell_trace,cell_struct(i).sen_stim_frames,opt.sta_pre_frames,opt.sta_post_frames,1:opt.sta_baseline_frames);
    if(~isempty(cell_struct(i).sen_stim_frames))
       cell_struct(i).sens_sta_amp = mean(cell_struct(i).sens_sta_traces(opt.sta_pre_frames:opt.sta_pre_frames+opt.sta_avg_frames));
       
        if  cell_struct(i).sens_sta_amp > opt.sta_thresh
            cell_struct(i).is_sensory = 1;
        end
       
    else
        cell_struct(i).sens_sta_amp = [];
    end
    % photostim sta
    if(this_num_trials>0)
        [~,~,~,~,~,cell_struct(i).photo_sta_traces,cell_struct(i).photo_sta_trace] = make_sta_from_traces(this_cell_trace,cell_struct(i).photostim_frames,opt.sta_pre_frames,opt.sta_post_frames,1:opt.sta_baseline_frames);
        cell_struct(i).photo_sta_amp = mean(cell_struct(i).photo_sta_trace(opt.sta_pre_frames:opt.sta_pre_frames+opt.sta_avg_frames));
        if  cell_struct(i).sta_amp > opt.sta_thresh
            cell_struct(i).is_photo = 1;
        end
   
        % plot traces
        hold on
        for t = 1:this_num_trials
            plot(cell_struct(i).photo_sta_traces(t,:),'color',opt.trial_color(t,:) ,'linewidth',1)
        end
        plot(cell_struct(i).photo_sta_trace,'color',[.5 .5 .5],'linewidth',1.5)
        plot([opt.sta_pre_frames opt.sta_pre_frames+photo_stim_duration],[0 0],'color','r','linewidth',2)

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
%%
figure('name','sens sta traces');
num_plot_cols = 6;
num_plot_rows = ceil(num_comp/num_plot_cols);
for i = 1:num_comp
    % plot traces
    subtightplot(num_plot_rows,num_plot_cols,i)
    hold on
    for t = 1:numel(cell_struct(i).sen_stim_frames)
        plot(cell_struct(i).sens_sta_traces(t,:),'color',opt.trial_color(t,:) ,'linewidth',1)
    end
    plot(cell_struct(i).sens_sta_trace,'color',[.5 .5 .5],'linewidth',1.5)
    ylimit = get(gca,'YLim');
    plot([opt.sta_pre_frames opt.sta_pre_frames+sen_stim_duration],[ylimit(1) ylimit(1)],'color',[.5 .5 .5],'linewidth',2) 
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

%%
figure('name','triggered sens sta traces');
num_plot_cols = 6;
num_plot_rows = ceil(num_comp/num_plot_cols);
for i = 1:num_comp
    % plot traces
    subtightplot(num_plot_rows,num_plot_cols,i)
    hold on
    if numel(cell_struct(i).trig_sen_stim_frames)>0
        for t = 1:numel(cell_struct(i).trig_sen_stim_frames)
            plot(cell_struct(i).trig_sen_sta_traces(t,:),'color',opt.trial_color(t,:) ,'linewidth',1)
        end
        plot(cell_struct(i).trig_sen_sta_trace,'color',[.5 .5 .5],'linewidth',1.5)
        ylimit = get(gca,'YLim');
        plot([opt.sta_pre_frames opt.sta_pre_frames+sen_stim_duration],[ylimit(1) ylimit(1)],'color',[.5 .5 .5],'linewidth',2) 
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
    if( cell_struct(i).is_sensory)
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