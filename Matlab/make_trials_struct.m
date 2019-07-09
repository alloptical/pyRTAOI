function [output_trials] = make_trials_struct(behavior_data)
% make trials.m struct from pybehavior output
fd_names = {'stim_type','correct','incorrect','fa','miss','firstresponse','cheated','auto_reward'};
num_fd = numel(fd_names);

output_trials = struct();
for i = 1:num_fd
    this_fd = fd_names{i};
    output_trials.(this_fd) = cellfun(@(x)extractfield(x,this_fd),behavior_data.results);
end

output_trials.('firstlick') = behavior_data.reaction_time';
%%
% figure
% test_mat = [output_trials.correct; output_trials.incorrect; output_trials.fa; output_trials.miss];
% imagesc(test_mat)
% unique(sum(test_mat,1))
%%
end

