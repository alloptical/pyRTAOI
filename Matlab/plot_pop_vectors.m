function plot_pop_vectors(pop_struct,plot_fields,num_states,opt,varargin)
% plot traces saved in get_pop_vectors
ylimit = [];
plot_ylabel = 'State probability';
IF_MEDIAN = 0;
plot_num_cols = 1;
IF_PLOT_RAW_ONLY = 0;
noise_thresh = [];
for v = 1:numel(varargin)
    if strcmpi(varargin{v},'ylimit')
        ylimit = varargin{v+1};
    end
    if strcmpi(varargin{v},'plot_ylabel')
        plot_ylabel = varargin{v+1};
    end
    if strcmpi(varargin{v},'IF_MEDIAN')
        IF_MEDIAN = varargin{v+1};
    end
    
    if strcmpi(varargin{v},'plot_num_cols')
        plot_num_cols = varargin{v+1};
    end
    
    if strcmpi(varargin{v},'IF_PLOT_RAW_ONLY')
        IF_PLOT_RAW_ONLY = varargin{v+1};
    end
    
    if strcmpi(varargin{v},'noise_thresh')
        noise_thresh = varargin{v+1};
    end
end

fds = plot_fields;
plot_num_rows = ceil(numel(fds)/plot_num_cols);
try
    try
        go_cue_frame = opt.go_cue_bin;
    catch
        try
            go_cue_frame = opt.gocue_bin;
        catch
            go_cue_frame = opt.sta_gocue_frame;
        end
        
    end
catch
    go_cue_frame = [];
end
try
    stim_on_frame = opt.sta_stimon_frame;
catch
    stim_on_frame = [];
end
state_colors = getOr(opt,'state_colors',brewermap(num_states,'Set1'));
trial_colors = getOr(opt,'trial_color',online_tex_init_color());
if num_states == 1
    state_colors = [.5 .5 .5];
end
if ~IF_PLOT_RAW_ONLY
    condi_colors = cell2mat(cellfun(@(f)trial_colors.(f),fds,'UniformOutput',false));
    if size(condi_colors,2)~=2
        condi_colors = cell2mat(cellfun(@(f)trial_colors.(f),fds,'UniformOutput',false)');
    end
    
    try
        if num_states>1
            figure('name','trial average')
            % imagesc
            for f = 1:numel(fds)
                this_F_avg = pop_struct.([fds{f} '_avg']);
                subplot(numel(fds),1,f)
                hold on
                imagesc(this_F_avg)
                plot([1 1].*go_cue_frame,ylim,'color','w','linewidth',2)
                ylim([0.5 0.5+num_states])
                xlabel('Time')
                ylabel('States')
                title(strrep(fds{f},'_',' '))
                
            end
        end
    end
    
    % shaded error bar
    figure('name','trial averge')
    for s = 1:num_states
        subplot(num_states,1,s)
        for f = 1:numel(fds)
            this_F = pop_struct.([fds{f}]);
            
            this_traces = this_F(:,:,s);
            if~isempty(this_traces)
                x_ticks =[0:1:size(this_traces,2)-1]./opt.Fs;
                hold on
                if ~IF_MEDIAN
                    % mean and sd
                    shadedErrorBar(x_ticks,mean(this_traces,1),...
                        std(this_traces,[],1),{'color',condi_colors(f,:),'linewidth',2},0.1);
                else
                    %         % median and quantile
                    shadedErrorBar(x_ticks,nanmedian(this_traces,1),[quantile(this_traces,0.75)-nanmedian(this_traces,1);...
                        nanmedian(this_traces,1)-quantile(this_traces,0.25)],{'color',condi_colors(f,:),'linewidth',2},0.5)
                    
                end
            end
        end
        if ~isempty(ylimit)
            ylim(ylimit)
        end
        plot([1 1].*go_cue_frame./opt.Fs,ylim,':','color','black','linewidth',2)
        plot(xlim,[0 0],'color','black')
        xlabel('Time')
        ylabel(plot_ylabel)
        title(['Component ' num2str(s) ' '])
        
    end
    
end
figure('name','raw projections','units','normalized','outerposition',[0 0 .6 1]); hold on
for f = 1:numel(fds)
    this_F = pop_struct.([fds{f}]);
    
    subplot(plot_num_rows,plot_num_cols,f)
    hold on
    % trace
    for s = 1:num_states
        plot(this_F(:,:,s)','color',state_colors(s,:));
    end
    if num_states == 1
        plot(mean(this_F(:,:,s),1),'color','black','linewidth',2);
    end
    
    % limits
    if ~isempty(ylimit)
        ylim(ylimit)
    end
    xlim([0, opt.trial_length-1])
    
    plot(xlim,[0 0],'color','r')
    axis square
    try
        plot([1 1].*go_cue_frame,ylim,':','color','black','linewidth',2)
        plot([1 1].*stim_on_frame,ylim,':','color','r','linewidth',2)
    end
    
    if~isempty(noise_thresh)
        plot(xlim,[1 1].*noise_thresh,'--','color','r')
        plot(xlim,-[1 1].*noise_thresh,'--','color','r')
        
    end
    
    xlabel('Time')
    ylabel(plot_ylabel)
    title([strrep(fds{f},'_',' ') ':  ' num2str(size(this_F,1)) ' trials'])
    
    
end
end

