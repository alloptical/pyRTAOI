function plot_pop_vectors(pop_struct,plot_fields,num_states,opt,varargin)
% plot traces saved in get_pop_vectors
ylimit = [];
plot_ylabel = 'State probability';
IF_MEDIAN = 0;
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
end

fds = plot_fields;
try
go_cue_frame = opt.go_cue_bin;
catch
    go_cue_frame = opt.gocue_bin;
end
state_colors = getOr(opt,'state_colors',brewermap(num_states,'Set1'));
trial_colors = getOr(opt,'trial_color',online_tex_init_color());
condi_colors = cell2mat(cellfun(@(f)trial_colors.(f),fds,'UniformOutput',false));
if size(condi_colors,2)~=2
    condi_colors = cell2mat(cellfun(@(f)trial_colors.(f),fds,'UniformOutput',false)');
end
if num_states == 1
    state_colors = [.5 .5 .5];
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
    if ~isempty(ylimit)
        ylim(ylimit)
    end
    plot([1 1].*go_cue_frame./opt.Fs,ylim,':','color','black','linewidth',2)
    xlabel('Time')
    ylabel(plot_ylabel)
    title(['Component ' num2str(s) ' '])
    
end


figure('name','trial raw')
for f = 1:numel(fds)
    this_F = pop_struct.([fds{f}]);
    subplot(numel(fds),1,f)
    hold on
    for s = 1:num_states
        plot(this_F(:,:,s)','color',state_colors(s,:));
    end
    plot([1 1].*go_cue_frame,ylim,':','color','black','linewidth',2)
    if ~isempty(ylimit)
        ylim(ylimit)
    end
    xlabel('Time')
    ylabel(plot_ylabel)
    xlim([0, size(this_F,2)])
    title([strrep(fds{f},'_',' ') ':  ' num2str(size(this_F,1)) ' trials'])
    if num_states == 1
        plot(mean(this_F(:,:,s),1),'color','black','linewidth',2); 
    end
end
end

