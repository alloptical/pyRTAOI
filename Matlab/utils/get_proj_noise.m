function [proj_sd,fds_of_interest] = get_proj_noise(choice_proj_struct,condition_stim_type,condition_stim_var,frames_of_interest)
% this is dummy
% get trajectory fluctualtions within-monitoring window 
fds_of_interest = {};
all_proj = [];
fd = 1;
for v = 1:length(condition_stim_type)
    this_stim = condition_stim_type(v);
    this_var = condition_stim_var(v);
    this_correct_fd =['stim_' num2str(this_stim) '_var_' num2str(this_var) '_correct'];
    this_incorrect_fd = ['stim_' num2str(this_stim) '_var_' num2str(this_var) '_incorrect'];
    this_projs = [choice_proj_struct.(this_correct_fd)(:,frames_of_interest);
        choice_proj_struct.(this_incorrect_fd)(:,frames_of_interest)];
    all_proj = [all_proj;this_projs];
    fds_of_interest{fd} = this_correct_fd;
    fds_of_interest{fd+1} = this_incorrect_fd;
    fd = fd+2;
end

proj_sd = mean(std(all_proj,[],2));
end

