function [baseline_cell_struct] = replot_get_cell_sta(baseline_cell_struct,this_sens_stim_frames,opt)
this_num_trials = numel(this_sens_stim_frames );
num_cells = size(baseline_cell_struct,2);
for i = 1:num_cells
    
    this_online_trace = baseline_cell_struct(i).filtC; % online trace med filtered
    
    baseline_cell_struct(i).sta_amp = 0;
    baseline_cell_struct(i).sta_traces = [];
    baseline_cell_struct(i).sta_trace = [];
    
    
    if(this_num_trials>0)
        [~,~,~,baseline_cell_struct(i).raw_sta_traces,~,~,baseline_cell_struct(i).raw_sta_trace] = make_sta_from_traces(this_online_trace,this_sens_stim_frames,opt.sta_pre_frames,opt.sta_post_frames,1:opt.sta_baseline_frames);
        % use df
%         [~,~,~,~,~,baseline_cell_struct(i).raw_sta_traces,baseline_cell_struct(i).raw_sta_trace] = make_sta_from_traces(this_online_trace,this_sens_stim_frames,opt.sta_pre_frames,opt.sta_post_frames,1:opt.sta_baseline_frames);
        
        if isfield(baseline_cell_struct,'deconvC')
            this_decov_trace = baseline_cell_struct(i).deconvC;
            [~,~,~,baseline_cell_struct(i).sta_traces,~,~,baseline_cell_struct(i).sta_trace] = make_sta_from_traces(this_decov_trace,this_sens_stim_frames,opt.sta_pre_frames,opt.sta_post_frames,1:opt.sta_baseline_frames);
%             [~,~,~,~,~,baseline_cell_struct(i).sta_traces,baseline_cell_struct(i).sta_trace] = make_sta_from_traces(this_decov_trace,this_sens_stim_frames,opt.sta_pre_frames,opt.sta_post_frames,1:opt.sta_baseline_frames);
            
        end
    end
    
end

end

