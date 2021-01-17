function [trans_struct] = plot_num_trans(vpath_struct,frames_of_interest,fds,trial_color)
for f = 1:numel(fds)
    fd = fds{f};
    this_traces = diff(vpath_struct.(fd),[],2);
    this_num_trials = size(this_traces,1);
    this_num_trans = nan(1,this_num_trials);
    for t = 1:this_num_trials
        this_num_trans(t) = length(find(this_traces(t,frames_of_interest)~=0)>0);
    end
    trans_struct.(fd) = this_num_trans;
end
condi_colors = cell2mat(cellfun(@(f)trial_color.(f)',fds,'UniformOutput',false))';
figure
scatter_cmp_conditions(trans_struct,[],...
        1,condi_colors,'connect_scatter',0,'BriefXlabel',0,'ShowMeanInXlabel',1,'add_jitter',1,'plot_stats',1);

end

