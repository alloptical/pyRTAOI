function [pop_struct,num_states] = get_pop_vectors(F,trial_length,trial_idx)
% F is a long population vector from concatinated trials (of equal length trial_length)
% trial_idx saves indices of different trials types in F  
pop_struct = struct();
fds = fields(trial_idx);
num_states = size(F,2);

for f = 1:numel(fds)
    trial_num = trial_idx.(fds{f});
    this_F = nan(length(trial_num),trial_length,num_states);
    for i = 1:length(trial_num)
        this_F(i,:,:) = F((trial_num(i)-1)*(trial_length)+[1:trial_length],:);
    end
    this_F_avg = squeeze(mean(this_F,1))';
    pop_struct.(fds{f}) = this_F;
%     pop_struct.([fds{f} '_avg']) = this_F_avg;
end

end

