function [cell_struct] = get_cell_active_auc(cell_struct,cmp_fds,opt,varargin)
% test if cells are turned off by stimulus
% compare frames with stim vs baseline frames
avg_frame_range = opt.sta_avg_frames;
baseline_avg_frames = opt.sta_baseline_frames+[1:30];
for v = 1:numel(varargin)
    if strcmpi(varargin{v},'baseline_avg_frames')
        baseline_avg_frames = varargin{v+1};
    end
end
for i = 1:size(cell_struct,2)
    for c = 1:numel(cmp_fds)
        
        all_stim = mean(cell_struct(i).(cmp_fds{c})(avg_frame_range,:),1);
        all_bs = mean(cell_struct(i).(cmp_fds{c})(baseline_avg_frames,:),1);
        
        labels = [ones(1,length(all_stim)),2.*ones(1,length(all_stim))]';
        scores = [ all_bs  all_stim]';
        [~,~,~, activeAUC] = perfcurve(labels,scores,2);
        cell_struct(i).([cmp_fds{c},'_activeAUC']) = activeAUC;
        if activeAUC >0.5
            cell_struct(i).([cmp_fds{c},'_offcell']) = 0;
        else
            cell_struct(i).([cmp_fds{c},'_offcell']) = 1;
            
        end
        
    end
    
end

