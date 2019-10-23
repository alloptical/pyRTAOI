function  plot_binary_decoder(stim_struct,stim_opt)
x_ticks = [1:stim_opt.trial_length]./stim_opt.Fs;
disc_frames_fd = getOr(stim_opt,'disc_frames_fd','shuf_disc_frames');
disc_frame_fd = getOr(stim_opt,'disc_frame_fd','shuf_disc_frame');
Nstd = getOr(stim_opt,'Nstd',2);
this_accuracy = stim_struct.classif_accuracy;
this_shuf_mean = stim_struct.shuf_classif_mean;
this_shuf_sd = stim_struct.shuf_classif_sd;
this_raw_trace = this_accuracy-this_shuf_mean;
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
plot([1,1].*stim_opt.stim_bin/stim_opt.Fs,ylim,':','color','r','linewidth',2)

plot([1,1].*dc_frame/stim_opt.Fs,ylim,'color','g','linewidth',1)

axis square
ylabel('Decoder accuracy (norm)')
xlabel('Time from trial start (sec)')



end

