function [] = plot_decoder_performance(deocd_struct,proj_struct,lick_struct,color,opt,varargin )
% adapted from 'plot_coding_direction' in BeastMatlabF\Texture dropbox
try
    trial_length = opt.trial_frames;
catch
    trial_length = opt.trial_length;
end
x_ticks = 1:trial_length;
zoomin_xlimit = [15,150];
Nstd = getOr(opt,'Nstd',1.5);

photo_color = [255 51 20]./255;
disc_frames_fd = getOr(opt,'disc_frames_fd','shuf_disc_frames');
disc_frame_fd = getOr(opt,'disc_frame_fd','shuf_disc_frame');
% plot_fds = {'st_correct_smooth_deconv_1','st_correct_smooth_deconv_2','st_incorrect_smooth_deconv_1','st_incorrect_smooth_deconv_2'};
fds =getOr(opt,'fd_names', fields(proj_struct));
correct_fd =fds(contains(fds, '_correct_'));
incorrect_fd =fds(contains(fds, '_incorrect_'));

plot_fds = {cell2mat(correct_fd(contains(correct_fd,'1'))),cell2mat(correct_fd(contains(correct_fd,'2')))...
    cell2mat(incorrect_fd(contains(incorrect_fd,'1'))),cell2mat(incorrect_fd(contains(incorrect_fd,'2')))};

this_colors = [color.correct_stim1; color.correct_stim2; color.incorrect_stim1; color.incorrect_stim2];


num_cols =5;
%% Raw traces with stim and lick
subplot(2,num_cols,1)
hold on
for f = [1,3]
    this_fd = plot_fds{f};
    this_trace =  proj_struct.(this_fd);
    plot( this_trace','color',this_colors(f,:))
    
    % lick time
    if ~isempty(lick_struct)
        this_lick_x= lick_struct.frame.(shorten_fd_name(this_fd));
        this_lick_y = nan(size(this_lick_x));
        for tr = 1:numel(this_lick_x)
            try
                this_lick_y(tr) = this_trace(tr,this_lick_x(tr));
            catch
                warning('lick frame out of trial frames!')
            end
        end
        scatter(this_lick_x,this_lick_y,30,'filled','MarkerEdgeColor',tint(this_colors(f,:),0.7), 'MarkerFaceColor',shade(this_colors(f,:),0.7))
    end
    % stim time
    try
    this_stim_x = deocd_struct.stim_frame.(this_fd);
    this_stim_y = nan(size(this_stim_x));
    for tr = 1:numel(this_stim_x)
        if ~isnan(this_stim_x(tr))
        try
            this_stim_y(tr) = this_trace(tr,this_stim_x(tr));
        catch
            warning('stim frame out of trial frames!')
        end
        end
    end
    scatter(this_stim_x,this_stim_y,40,'diamond','filled','MarkerEdgeColor',tint(photo_color,0.7), 'MarkerFaceColor',photo_color)
    catch
        warning('no photostim found')
    end
    
end
try
    photo_frame = deocd_struct.('stim1').photo_enable_frame ;
    hit_rate = deocd_struct.stim1.hit_rate;
    fa_rate = deocd_struct.stim1.fa_rate;
    fscore =  deocd_struct.stim1.fscore;
    
    plot(deocd_struct.thresh_framewise.('stim1'),'color','black','linewidth',2)
    plot([1,1].*photo_frame(1),ylim,'color',photo_color,'linewidth',2)
    plot([1,1].*photo_frame(end),ylim,':','color',photo_color,'linewidth',2)
   
    % mark discriminating frames
    % disc_frames = cd_struct.stim1.disc_frames;
    % scatter(disc_frames,zeros(size(disc_frames)),30,'square','filled','MarkerEdgeColor','none', 'MarkerFaceColor','black')
catch
    hit_rate = 0;
    fa_rate = 1;
    fscore =  0;
end

plot([1,1].*opt.gocue_frame_adj,ylim,':','color',[.5 .5 .5],'linewidth',2)
xlim(zoomin_xlimit)
axis square  

xlabel('Frames')
ylabel('Projection')
title({['(before lick) hit: ' num2str(hit_rate*100,'%10.2f') '% ,fa:' num2str(fa_rate*100,'%10.2f') '%']; ['F1 score:'  num2str(fscore,'%10.2f')]})

subplot(2,num_cols,num_cols+1)
hold on
for f = [2,4]
    this_fd = plot_fds{f};
    this_trace =  proj_struct.(this_fd);
    plot( this_trace','color',this_colors(f,:))
    
    % lick time
    if ~isempty(lick_struct)
        this_lick_x= lick_struct.frame.(shorten_fd_name(this_fd));
        this_lick_y = nan(size(this_lick_x));
        for tr = 1:numel(this_lick_x)
            try
                this_lick_y(tr) = this_trace(tr,this_lick_x(tr));
            catch
                warning('lick frame out of trial frames!')
            end
        end
        scatter(this_lick_x,this_lick_y,30,'filled','MarkerEdgeColor',tint(this_colors(f,:),0.7), 'MarkerFaceColor',shade(this_colors(f,:),0.7))
    end
    % stim time
    try
    this_stim_x= deocd_struct.stim_frame.(this_fd);
    this_stim_y = nan(size(this_stim_x));
    for tr = 1:numel(this_stim_x)
        if ~isnan(this_stim_x(tr))
            try
                this_stim_y(tr) = this_trace(tr,this_stim_x(tr));
            catch
                warning('stim frame out of trial frames!')
            end
        end
    end
    scatter(this_stim_x,this_stim_y,40,'diamond','filled','MarkerEdgeColor',tint(photo_color,0.7), 'MarkerFaceColor',photo_color)
    catch
        warning('no photostim found')
    end
    
end
try
    photo_frame =  deocd_struct.('stim2').photo_enable_frame ;
    hit_rate = deocd_struct.stim2.hit_rate;
    fa_rate = deocd_struct.stim2.fa_rate;
    fscore =  deocd_struct.stim2.fscore;
    
    plot(deocd_struct.thresh_framewise.('stim2'),'color','black','linewidth',2)
    plot([1,1].*photo_frame(1),ylim,'color',photo_color,'linewidth',2)
    plot([1,1].*photo_frame(end),ylim,':','color',photo_color,'linewidth',2)
    
    % mark discriminating frames
    % disc_frames = cd_struct.stim1.disc_frames;
    % scatter(disc_frames,zeros(size(disc_frames)),30,'square','filled','MarkerEdgeColor','none', 'MarkerFaceColor','black')
catch
    hit_rate = 0;
    fa_rate = 1;
    fscore =  0;
end

plot([1,1].*opt.gocue_frame_adj,ylim,':','color',[.5 .5 .5],'linewidth',2)
xlim(zoomin_xlimit)
axis square  

xlabel('Frames')
ylabel('Projection')
title({['hit: ' num2str(hit_rate*100,'%10.2f') '% ,fa:' num2str(fa_rate*100,'%10.2f') '%']; ['F1 score:'  num2str(fscore,'%10.2f')]})


%% SHADED ERROR BAR PLOTS
subplot(2,num_cols,2)
hold on
for f = [1,3]
    this_fd =plot_fds{f};
    temp = strsplit(this_fd,'_');
    this_traces =  proj_struct.(this_fd);
    this_thresh = deocd_struct.thresh_framewise.(['stim' temp{end}]);
    
    for tr = 1:size(this_traces,1)
        this_traces(tr,:) = this_traces(tr,:) - this_thresh;
    end
    
    
%     shadedErrorBar(x_ticks,mean(this_traces,1),...
%         std(this_traces,[],1),{'color',color.(this_fd),'linewidth',2},0.1)
    shadedErrorBar(x_ticks,nanmedian(this_traces,1),[quantile(this_traces,0.75)-nanmedian(this_traces,1);...
        nanmedian(this_traces,1)-quantile(this_traces,0.25)],{'color',this_colors(f,:),'linewidth',2},0.5);
    
end
% plot(cd_struct.cd_thresh_log.('stim1'),'--','color','black','linewidth',2)
plot(xlim,[0 0],'color','black')
plot([1,1].*opt.gocue_frame_adj,ylim,':','color',[.5 .5 .5],'linewidth',2)

% mark discriminating frames
disc_frames = deocd_struct.stim1.(disc_frames_fd);
scatter(disc_frames,zeros(size(disc_frames)),30,'square','filled','MarkerEdgeColor','none', 'MarkerFaceColor','black')

xlim([0,x_ticks(end)])
axis square  
text(0.1,1,['Num correct trials: ' num2str(size(proj_struct.(plot_fds{1}),1))],...
    'units','normalized', 'horizontalalignment','left', 'color',this_colors(1,:))
text(0.1,0.9,['Num incorrect trials: ' num2str(size(proj_struct.(plot_fds{3}),1))],...
    'units','normalized', 'horizontalalignment','left', 'color',this_colors(3,:))
xlabel('Frames')
ylabel('Projection')

subplot(2,num_cols,num_cols+2)
hold on
for f = [2,4]
    this_fd =plot_fds{f};
    temp = strsplit(this_fd,'_');
    this_traces =  proj_struct.(this_fd);%- cd_struct.cd_thresh_log.(['stim' temp{end}]);
    this_thresh = deocd_struct.thresh_framewise.(['stim' temp{end}]);
    
    for tr = 1:size(this_traces,1)
        this_traces(tr,:) = this_traces(tr,:) - this_thresh;
    end
%     shadedErrorBar(x_ticks,mean(this_traces,1),...
%         std(this_traces,[],1),{'color',color.(this_fd),'linewidth',2},0.1)
    shadedErrorBar(x_ticks,median(this_traces,1),[quantile(this_traces,0.75)-median(this_traces,1);...
        median(this_traces,1)-quantile(this_traces,0.25)],{'color',this_colors(f,:),'linewidth',2},0.5);
%    plot(xlim,[1,1].*cd_struct.cd_thresh_log.('stim1'),'--','color','black','linewidth',2)
 
end

plot(xlim,[0 0],'color','black')
plot([1,1].*opt.gocue_frame_adj,ylim,':','color',[.5 .5 .5],'linewidth',2)

% mark discriminating frames
disc_frames = deocd_struct.stim2.(disc_frames_fd);
scatter(disc_frames,zeros(size(disc_frames)),30,'square','filled','MarkerEdgeColor','none', 'MarkerFaceColor','black')

xlim([0,x_ticks(end)])
axis square  


text(0.1,1,['Num correct trials: ' num2str(size(proj_struct.(plot_fds{2}),1))],...
    'units','normalized', 'horizontalalignment','left', 'color',this_colors(2,:))
text(0.1,0.9,['Num incorrect trials: ' num2str(size(proj_struct.(plot_fds{4}),1))],...
    'units','normalized', 'horizontalalignment','left', 'color',this_colors(4,:))

xlabel('Frames')
ylabel('Projection')
%% classification accuracy
min_frames = 15;
subplot(2,num_cols,3)
stim_fd = 'stim1';
hold on
this_accuracy = deocd_struct.(stim_fd).classif_accuracy;
this_shuf_mean = deocd_struct.(stim_fd).shuf_classif_mean;
this_shuf_sd = deocd_struct.(stim_fd).shuf_classif_sd;
this_raw_trace = this_accuracy-this_shuf_mean;
shadedErrorBar(x_ticks,zeros(1,opt.trial_frames),...
    Nstd.*this_shuf_sd,{'color',[.5 .5 .5],'linewidth',2},0.5);
plot(this_raw_trace,'color','black')
this_dc_frame =  deocd_struct.stim1.(disc_frame_fd);
dc_frame_stim1 = this_dc_frame;
ylim([-0.5 0.5])
plot([1,1].*opt.gocue_frame_adj,ylim,':','color',[.5 .5 .5],'linewidth',2)


try
    plot([1 1].*this_dc_frame,ylim,':','color','r','linewidth',1.5)
    text(1,0.1,['t = ' num2str(this_dc_frame/opt.frame_rate,'%10.2f'), 's'], 'units','normalized', 'horizontalalignment','right', 'color',[1 0 0])
end
ylabel('Norm. accuracy')
xlabel('Frames')
axis square

subplot(2,num_cols,num_cols+3)
stim_fd = 'stim2';
hold on
this_accuracy = deocd_struct.(stim_fd).classif_accuracy;
this_shuf_mean = deocd_struct.(stim_fd).shuf_classif_mean;
this_shuf_sd = deocd_struct.(stim_fd).shuf_classif_sd;
this_raw_trace = this_accuracy-this_shuf_mean;
shadedErrorBar(x_ticks,zeros(1,opt.trial_frames),...
    Nstd.*this_shuf_sd,{'color',[.5 .5 .5],'linewidth',2},0.5);
plot(this_raw_trace,'color','black')
this_dc_frame =  deocd_struct.stim2.(disc_frame_fd);
dc_frame_stim2 = this_dc_frame;

ylim([-0.5 0.5])
plot([1,1].*opt.gocue_frame_adj,ylim,':','color',[.5 .5 .5],'linewidth',2)


try
    this_dc_frame = this_dc_frame(1);
    plot([1 1].*this_dc_frame,ylim,':','color','r','linewidth',1.5)
    text(1,0.1,['t = ' num2str(this_dc_frame/opt.frame_rate,'%10.2f'), 's'], 'units','normalized', 'horizontalalignment','right', 'color',[1 0 0])
end
ylabel('Norm. accuracy')
xlabel('Frames')
axis square

%% framewise hit and fa rate
subplot(2,num_cols,4)
stim_fd = 'stim1';

hold on
hit_rate =  deocd_struct.(stim_fd).framewise_hr;
fa_rate = deocd_struct.(stim_fd).framewise_fa;
plot(x_ticks,hit_rate,'color',[153 255 153]./255)
plot(x_ticks,fa_rate,'color',[255 153 153]./255)
try
    this_dc_frame = dc_frame_stim1(1);
    plot([1 1].*this_dc_frame,[0 1],':','color','r','linewidth',1.5)
    text(1,0.1,['t = ' num2str(this_dc_frame/opt.frame_rate,'%10.2f'), 's'], 'units','normalized', 'horizontalalignment','right', 'color',[1 0 0])
end
ylim([0 1])
xlabel('Frames')
ylabel('Hit and FA rate')
axis square



subplot(2,num_cols,num_cols+4)
stim_fd = 'stim2';

hold on
hit_rate =  deocd_struct.(stim_fd).framewise_hr;
fa_rate = deocd_struct.(stim_fd).framewise_fa;
plot(x_ticks,hit_rate,'color',[153 255 153]./255)
plot(x_ticks,fa_rate,'color',[255 153 153]./255)
try
    this_dc_frame = dc_frame_stim2(1);
    plot([1 1].*this_dc_frame,[0 1],':','color','r','linewidth',1.5)
    text(1,0.1,['t = ' num2str(this_dc_frame/opt.frame_rate,'%10.2f'), 's'], 'units','normalized', 'horizontalalignment','right', 'color',[1 0 0])
end
ylim([0 1])
xlabel('Frames')
ylabel('Hit and FA rate')
axis square



%% Detection time and lick time (contineu check photostim)
subplot(2,num_cols,5)
try
    values = struct();
    values.stim = deocd_struct.stim1.this_det_times;
    values.lick = deocd_struct.stim1.this_lick_time;
    scatter_cmp_conditions(values,[],...
        1,[photo_color;0 0 0],'connect_scatter',1,'BriefXlabel',1,'ShowMeanInXlabel',1,'add_jitter',0);
    ylabel('Time from trial start (sec)')
catch
    text(0,1,['Stim1 outcome prediction error'],...
        'units','normalized', 'horizontalalignment','left', 'color','r')
    axis off
end
axis square


subplot(2,num_cols,5+num_cols)
try
values = struct();
values.stim = deocd_struct.stim2.this_det_times;
values.lick = deocd_struct.stim2.this_lick_time;
scatter_cmp_conditions(values,[],...
    1,[photo_color;0 0 0],'connect_scatter',1,'BriefXlabel',1,'ShowMeanInXlabel',1,'add_jitter',0);
ylabel('Time from trial start (sec)')
axis square  
catch
    text(0.1,1,['Stim2 outcome prediction error'],...
        'units','normalized', 'horizontalalignment','left', 'color','r')
    axis off
end
axis square

end

