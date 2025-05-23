function [state_idx_struct,corr_struct,FLAG_MATCHED_STATES] = get_state_idx(hmm_struct,epoch_window,ref_fields,num_states)
% see if there is a state significantly more correlated with
% texture-related window 
num_bins = size(hmm_struct.(ref_fields{1}),2);
ref_trace = zeros(1,num_bins);
ref_trace(epoch_window) = 1; % reference trace where decision bins are ones
corr_struct = struct();
FLAG_MATCHED_STATES = zeros(1,numel(ref_fields));
for f = 1:numel(ref_fields)
    this_fd = ref_fields{f};
    state_idx_struct.(this_fd) = 0;
    this_traces = hmm_struct.(this_fd);
    this_num_trials = size(this_traces,1);
    corrs = nan(this_num_trials,num_states);

    for s = 1:num_states
        for t = 1:this_num_trials % loop through trials
%             corrs(t,s)  = xcorr(ref_trace',this_traces(t,:,s)',0,'coeff'  );% correlation with state trial average  - dosent make sense  
            corrs(t,s)  = sum(ref_trace.*this_traces(t,:,s))/length(epoch_window); % probability

        end
    end
    [p,tbl,stats] = anova1(corrs,[],'off'); % anova test (colums are groups)
	[comparison,means,h,gnames] = multcompare(stats,'display','off');  
    if   all(comparison(:,6)<0.01)
        FLAG_MATCHED_STATES(f) = find(means(:,1)==max(means(:,1)));
        state_idx_struct.(this_fd) = FLAG_MATCHED_STATES(f);
    else
        for s = 1:num_states
            this_p = comparison(find(comparison(:,1)==s|comparison(:,2)==s),end);
%             if p<0.05 && all(this_p<0.01/size(comparison,1)) % too strict
            if p<0.05 && numel(find(this_p<0.05/size(comparison,1)))>round(numel(this_p)/2) % significant among more than half of comprisons
                FLAG_MATCHED_STATES(f) = s;
                state_idx_struct.(this_fd) = s;
            end
        end
    end
    
    % save the state with max mean if didnt pass stats test
    if  FLAG_MATCHED_STATES(f)==0
         [~,state_idx_struct.(this_fd)] = max(stats.means);
    end
    


    corr_struct.(this_fd).data = corrs;
    corr_struct.(this_fd).stats = stats;
end

if all(FLAG_MATCHED_STATES>0)
    disp('Found matched states with trial structure')
    disp(state_idx_struct)
else
    warning('Trial-related states not found! ')
end
end

