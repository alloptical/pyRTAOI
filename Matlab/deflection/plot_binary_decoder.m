function  plot_binary_decoder(stim_struct,stim_opt,varargin)
x_ticks = [1:stim_opt.trial_length]./stim_opt.Fs;
x_label = 'Time from trial start (sec)';
IF_ALIGNGOCUE = false;
for v = 1:numel(varargin)
    if strcmpi(varargin{v},'IF_ALIGNGOCUE')
        IF_ALIGNGOCUE = varargin{v+1};
    end
    
    if strcmpi(varargin{v},'x_label')
        x_label = varargin{v+1};
    end
    
end
disc_frames_fd = getOr(stim_opt,'disc_frames_fd','shuf_disc_frames');
disc_frame_fd = getOr(stim_opt,'disc_frame_fd','shuf_disc_frame');
Nstd = getOr(stim_opt,'Nstd',2);
this_accuracy = stim_struct.classif_accuracy;
this_shuf_mean = stim_struct.shuf_classif_mean;
this_shuf_sd = stim_struct.shuf_classif_sd;
this_raw_trace = this_accuracy-this_shuf_mean;

subplot(1,2,1); hold on
shadedErrorBar(x_ticks,zeros(1,length(x_ticks)),...
    Nstd.*this_shuf_sd,{'color',[.5 .5 .5],'linewidth',2},0.5);
plot(x_ticks,this_raw_trace,'color','black')
this_dc_frame =  stim_struct.(disc_frame_fd);
this_dc_frames =  stim_struct.(disc_frames_fd);

dc_frame = this_dc_frame;
scatter(this_dc_frames/stim_opt.Fs,zeros(size(this_dc_frames)),30,'square','filled','MarkerEdgeColor','none', 'MarkerFaceColor','black')

ylim([-0.5 0.5])
xlim([x_ticks(1),x_ticks(end)])
plot([1,1].*stim_opt.gocue_bin/stim_opt.Fs,ylim,':','color',[.5 .5 .5],'linewidth',2)
% try
% plot([1,1].*stim_opt.stim_bin/stim_opt.Fs,ylim,':','color','r','linewidth',2)
% end
plot([1,1].*dc_frame/stim_opt.Fs,ylim,'color','g','linewidth',1)

axis square
ylabel('Decoder accuracy (norm)')
xlabel(x_label)
plot_ticks = xticks;
if IF_ALIGNGOCUE
    xticklabels(plot_ticks-stim_opt.gocue_bin/stim_opt.Fs)
% set(gca, 'XTickLabel', x_tick_label)
end

subplot(1,2,2); hold on
plot(x_ticks,stim_struct.framewise_hr,'color','black')
plot(x_ticks,stim_struct.framewise_fa,'color',[.5 .5 .5])
axis square
ylim([0 1])
xlim([x_ticks(1),x_ticks(end)])
plot([1,1].*stim_opt.gocue_bin/stim_opt.Fs,ylim,':','color',[.5 .5 .5],'linewidth',2)
% try
% plot([1,1].*stim_opt.stim_bin/stim_opt.Fs,ylim,':','color','r','linewidth',2)
% end
try
plot([1,1].*stim_struct.disc_frame/stim_opt.Fs,ylim,'color','g','linewidth',1)
end
legend('Hit','FA')
plot_ticks = xticks;
if IF_ALIGNGOCUE
    xticklabels(plot_ticks-stim_opt.gocue_bin./stim_opt.Fs)
end
end

