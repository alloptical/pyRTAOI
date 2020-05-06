function [] = plot_trial_avg_imagesc(cell_struct,trial_indices,photo_ensembles,plot_cell_idx,plot_photo_types,plot_condition_types,opt,varargin)
sta_trace_fd = 'sta_traces';
sort_cell = true;
for v = 1:numel(varargin)
    if strcmpi(varargin{v},'sta_trace_fd')
        sta_trace_fd = varargin{v+1};
    end
    if strcmpi(varargin{v},'sort_cell')
        sort_cell = varargin{v+1};
    end
end

% compare three types of trial avg
num_photo_types = numel(plot_photo_types);
num_compares = nchoosek(num_photo_types,2);
cmp_pairs = combnk(1:num_photo_types,2);

num_plot_cols = num_photo_types+num_compares;
num_plot_rows = 2;
all_trial_types = fields(trial_indices);

fun = @(x)text(0,x,['\rightarrow'], 'horizontalalignment','right', 'color','red');
xaxisvalues = [0:30:opt.sta_pre_frames+opt.sta_post_frames];
figure('name',['cell avg traces'],'units','normalized','outerposition',[0 0 1 1]);
plot_count = 1;

for p = 1:numel(plot_condition_types)
    photo_idx = plot_condition_types(p);
    all_traces = {};
    all_trials = {};
    for photo_type = 1:numel(plot_photo_types)
        this_fds = all_trial_types(contains(all_trial_types,['stim_' num2str(photo_idx) plot_photo_types{photo_type}]));
        this_trials = unique(cell2mat(cellfun(@(x)trial_indices.(x)(:),this_fds,'UniformOutput',false)));
        ax = subplot(num_plot_rows,num_plot_cols,plot_count);hold on
        
        if ~isempty(this_trials)
            this_traces = cell2mat(arrayfun(@(x)mean(x.(sta_trace_fd)(this_trials,:),1)',cell_struct(plot_cell_idx),'UniformOutput',false));
            
            if photo_type ==1
                if sort_cell
                % rank cell by sta amp and set colorlut limit
                this_amp = mean(this_traces(opt.sta_avg_frames,:),1);
                [~,this_plot_cell_idx] = sort(this_amp,'descend');
                ylabel('Cell idx (sorted)')
                else
                    this_plot_cell_idx = plot_cell_idx;
                    ylabel('Cell idx')

                end
            end
            
            if plot_count ==1
                zlimit =[min(this_traces(:)) max(this_traces(:))];
            end
            
            if photo_type == 1
                try
                    [~,mark_target_idx] = intersect(this_plot_cell_idx,photo_ensembles{photo_idx});
                catch
                    mark_target_idx = [];
                end
            end
            
            this_traces = this_traces(:,this_plot_cell_idx)';
            imagesc(this_traces)
            colormap(ax,b2r(zlimit(1) ,zlimit(2)))
            
            % mark target idx
            arrayfun(@(x)fun(x),mark_target_idx)
            
            % mark go-cue
            plot([1,1].*opt.sta_gocue_frame, ylim,'color','black','linestyle',':')
            axis ij
            box off
            ylim([0 size(this_traces,2)])
            xlim([0 opt.sta_pre_frames+opt.sta_post_frames+1])
            xticks(xaxisvalues)
            xticklabels(arrayfun(@(x){num2str(x)},(xaxisvalues-opt.sta_trialon_frame)./opt.frame_rate))
        else
            this_traces = [];
        end
        plot_count = plot_count+1;
%         colorbar('SouthOutside')
        
        title({['Tex' num2str(photo_idx) ' '  strrep(plot_photo_types{photo_type},'_',' ') ];[ num2str(numel(this_trials)) ' trials']})
        
        all_traces{photo_type} = this_traces;
        all_trials{photo_type} = this_trials;
    end
%     colorbar('EastOutside')
   
    % differences
    for c = 1:num_compares
        c1 = cmp_pairs(c,1);
        c2 = cmp_pairs(c,2);
        ax = subplot(num_plot_rows,num_plot_cols,plot_count); hold on

        if (~isempty(all_traces{c1}))&&(~isempty(all_traces{c2}))
            
            diff_traces = all_traces{c1}-all_traces{c2};
            if c == 1&& p ==1
                diff_zlimit =[min(diff_traces(:)) max(diff_traces(:))];
            end
            imagesc(diff_traces);
            colormap(ax,b2r(diff_zlimit(1) ,diff_zlimit(2)))
            arrayfun(@(x)fun(x),mark_target_idx)
            plot([1,1].*opt.sta_gocue_frame, ylim,'color','black','linestyle',':')
            axis ij
            box off
            ylim([0 size(diff_traces,2)])
            
            xlim([0 opt.sta_pre_frames+opt.sta_post_frames+1])
            xticks(xaxisvalues)
            xticklabels(arrayfun(@(x){num2str(x)},(xaxisvalues-opt.sta_trialon_frame)./opt.frame_rate))
        end
        plot_count = plot_count+1;
        title({[strrep(plot_photo_types{c1},'_',' ')];[ '-' strrep(plot_photo_types{c2},'_',' ')]})
    end
    xlabel('Time from trial start (sec)')
    %     this_title = ['TrialAvgStim' num2str(photo_idx) '_STATrace_' strrep(caiman_file,'.mat','')];
    %     suptitle(strrep(this_title,'_',' '))
    %     export_fig([fig_save_path filesep this_title '.png'])
    
end

end

