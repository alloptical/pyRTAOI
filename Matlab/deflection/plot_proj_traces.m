function  plot_proj_traces(pop_struct,fds,opt,varargin)
x_ticks = [1:opt.trial_length]./opt.Fs;
ylimit = [];
plot_ylabel = 'Projection';
for v = 1:numel(varargin)
    if strcmpi(varargin{v},'ylimit')
        ylimit = varargin{v+1};
    end
    if strcmpi(varargin{v},'plot_ylabel')
        plot_ylabel = varargin{v+1};
    end
end
trial_colors = getOr(opt,'trial_color',deflect_init_color());
condi_colors = cell2mat(cellfun(@(f)trial_colors.(f),fds,'UniformOutput',false));
if size(condi_colors,2)~=2
    condi_colors = cell2mat(cellfun(@(f)trial_colors.(f),fds,'UniformOutput',false)');
end
try
    go_cue_frame = opt.go_cue_bin;
catch
    go_cue_frame = opt.gocue_bin;
end
for f = 1:numel(fds)
    this_traces = pop_struct.([fds{f}]);
    hold on
    % mean and sd
    shadedErrorBar(x_ticks,mean(this_traces,1),...
        std(this_traces,[],1),{'color',condi_colors(f,:),'linewidth',2},0.1);
    
    %         % median and quantile
    %         shadedErrorBar(x_ticks,nanmedian(this_traces,1),[quantile(this_traces,0.75)-nanmedian(this_traces,1);...
    %             nanmedian(this_traces,1)-quantile(this_traces,0.25)],{'color',condi_colors(f,:),'linewidth',2},0.5)
end
if ~isempty(ylimit)
    ylim(ylimit)
end
plot([1 1].*go_cue_frame./opt.Fs,ylim,':','color','black','linewidth',2)
xlabel('Time')
ylabel(plot_ylabel)
    
xlim([x_ticks(1),x_ticks(end)])
plot([1,1].*opt.gocue_bin/opt.Fs,ylim,':','color',[.5 .5 .5],'linewidth',2)
plot([1,1].*opt.stim_bin/opt.Fs,ylim,':','color','r','linewidth',2)


axis square
ylabel('Projection')
xlabel('Time from trial start (sec)')



end

