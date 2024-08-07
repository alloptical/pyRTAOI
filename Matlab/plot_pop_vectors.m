function plot_pop_vectors(pop_struct,plot_fields,num_states,opt,varargin)
% plot traces saved in get_pop_vectors
ylimit = [];
xlimit = [0, opt.trial_length-1];
plot_ylabel = 'State probability';
IF_MEDIAN = 0;
plot_num_cols = 1;
IF_PLOT_RAW_ONLY = 0;
IF_PLOT_AVG_ONLY = 0;
IF_SAVE_PLOT = 0;
IF_NORMALISE = false;
noise_thresh = [];
plot_area = false; % shaded area under curve, for hmm states
sup_title = [];
func = @(x)1./(1+exp(-x))-0.5;
for v = 1:numel(varargin)
    if strcmpi(varargin{v},'ylimit')
        ylimit = varargin{v+1};
    end
    if strcmpi(varargin{v},'xlimit')
        xlimit = varargin{v+1};
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
    
    if strcmpi(varargin{v},'plot_area')
        plot_area = varargin{v+1};
    end
    
    if strcmpi(varargin{v},'sup_title')
        sup_title = varargin{v+1};
    end
    
    if strcmpi(varargin{v},'IF_PLOT_AVG_ONLY')
        IF_PLOT_AVG_ONLY = varargin{v+1};
    end
    
    if strcmpi(varargin{v},'IF_SAVE_PLOT')
        IF_SAVE_PLOT = varargin{v+1};
    end
    
    if strcmpi(varargin{v},'IF_NORMALISE')
        IF_NORMALISE = varargin{v+1}; %logistic regression, normalise to max probability==1
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
try
    opt.Fs = opt.frame_rate;
end
state_colors = getOr(opt,'state_colors',brewermap(num_states,'Set1'));
trial_colors = getOr(opt,'trial_color',online_tex_init_color());
if num_states == 1
    state_colors = [.5 .5 .5];
end
if ~IF_PLOT_RAW_ONLY
    try
    condi_colors = cell2mat(cellfun(@(f)trial_colors.(strrep(f,'_smooth_deconv','_stim')),fds,'UniformOutput',false));
    catch
            condi_colors = cell2mat(cellfun(@(f)trial_colors.(f),fds,'UniformOutput',false));

    end
    if size(condi_colors,2)~=3
        condi_colors = cell2mat(cellfun(@(f)trial_colors.(strrep(f,'_smooth_deconv','_stim')),fds,'UniformOutput',false)');
    end
    

    % shaded error bar
    figure('name','trial averge','units','normalized','outerposition',[0 0 .4,.7])
    plot_count = 1;
    for s = 1:num_states
        subplot(num_states,2,plot_count); hold on
        for f = 1:numel(fds)
            try
            this_F = pop_struct.([fds{f}]);
            if IF_NORMALISE
                this_F = func(this_F);
            end
            this_traces = this_F(:,:,s);

            plot(this_traces','color',condi_colors(f,:)');
            catch
                continue
            end
        end
        if ~isempty(ylimit)
            ylim(ylimit)
        end

        xlim(xlimit)
        plot([1 1].*go_cue_frame,ylim,':','color','black','linewidth',2)
        plot(xlim,[0 0],'color','black')
        xlabel('Time')
        ylabel(plot_ylabel)
        if num_states>1
        title(['Component ' num2str(s) ' '])
        end
        plot_count = plot_count+1;
        axis square
        
        subplot(num_states,2,plot_count)
        for f = 1:numel(fds)
            try
            this_F = pop_struct.([fds{f}]);
            if IF_NORMALISE
                this_F = func(this_F);
            end
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
            catch
                continue
            end
        end
        if ~isempty(ylimit)
            ylim(ylimit)
        end
        plot([1 1].*go_cue_frame./opt.Fs,ylim,':','color','black','linewidth',2)
        plot(xlim,[0 0],'color','black')
        xlabel('Time')
        ylabel(plot_ylabel)
        if num_states>1
        title(['Component ' num2str(s) ' '])
        end
        plot_count = plot_count+1;
        axis square

        
        
    end
    if ~isempty(sup_title)
        suptitle(strrep(sup_title,'_',''))
    end
    if IF_SAVE_PLOT
       export_fig([opt.save_path filesep sup_title '_AVG.png'])
    end

    
end
if ~IF_PLOT_AVG_ONLY
    figure('name','raw projections','units','normalized','outerposition',[0 0 1 1]); hold on
    for f = 1:numel(fds)
        try
        this_F = pop_struct.([fds{f}]);
        if IF_NORMALISE
            this_F = func(this_F);
        end
        catch
            disp([fds{f} 'error, skipped'])
            continue
        end
        
        subplot(plot_num_rows,plot_num_cols,f)
        hold on
        % trace
        for s = 1:num_states
            if plot_area
                for t = 1:size(this_F,1)
                    area(this_F(t,:,s)','FaceColor',tint(state_colors(s,:),0.5),'FaceAlpha',0.1,'EdgeColor','none');
                end
            end
            plot(this_F(:,:,s)','color',state_colors(s,:));
            
        end
        if num_states == 1
            plot(mean(this_F(:,:,s),1),'color','black','linewidth',2);
        end
        
        % limits
        if ~isempty(ylimit)
            ylim(ylimit)
        end
          xlim(xlimit)

        
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
        title([strrep(strrep(fds{f},'deconv',''),'_',' ') ': ' num2str(size(this_F,1))])
        
    end
    if ~isempty(sup_title)
        suptitle(strrep(sup_title,'_',''))
    end
    
    if IF_SAVE_PLOT
       export_fig([opt.save_path filesep sup_title '_RAW.png'])
    end
end
end

