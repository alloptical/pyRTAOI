function [output_trials,odd_trial_idx] = make_trials_struct(behavior_data)
% make trials.m struct from pybehavior output
fd_names = {'stim_type','correct','incorrect','miss','firstresponse','cheated','auto_reward','stim_var','fa'};
num_fd = numel(fd_names);
odd_trial_idx = [];
output_trials = struct();
tot_num_trials = numel(behavior_data.results);
for i = 1:num_fd
    this_fd = fd_names{i};
    if isfield(behavior_data.results{1},this_fd)
        try
        this_data = cellfun(@(x)extractfield(x,this_fd),behavior_data.results,'UniformOutput',false);
        catch
                    this_data = cellfun(@(x)extractfield(x,this_fd),behavior_data.results(1:tot_num_trials-1),'UniformOutput',false);
        end
        if iscell(this_data)
            checktype = cellfun(@(x)class(x),this_data,'UniformOutput',false);
            [types,~,whichcell]=unique(checktype);
            if contains(types,'cell')
                this_data= cellfun(@(x)cell2mat(x),this_data,'UniformOutput',false);
            end
            [this_type_idx]= mode(whichcell);
            this_data = cell2mat(this_data(whichcell== this_type_idx));

            this_odd_trial = find(whichcell~= this_type_idx);
            odd_trial_idx = [odd_trial_idx,this_odd_trial];
            if numel(types)>1
                warning(['some trial is odd! trial ' num2str(this_odd_trial) ' was excluded in ' this_fd])
            end
            
        end
        output_trials.(this_fd) = this_data;
    else
        warning(['field not found in behavior: ' this_fd ', skipped'])
    end
end
output_trials.('firstlick') = behavior_data.reaction_time(1:tot_num_trials)';
% exclude odd trials
odd_trial_idx = unique(odd_trial_idx);
normal_trial_idx = setdiff(1:tot_num_trials,odd_trial_idx);
if ~isempty(odd_trial_idx)
    disp(['odd trials indices: ' num2str(odd_trial_idx) ', excluded']) % will get such trials if pybehav was forced to abort
    output_trials = structfun(@(x)x(normal_trial_idx),output_trials,'UniformOutput',false);
end
%%
% figure
% test_mat = [output_trials.correct; output_trials.incorrect; output_trials.fa; output_trials.miss];
% imagesc(test_mat)
% unique(sum(test_mat,1))
%%
end

