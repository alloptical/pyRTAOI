function [pop_struct,num_states] = get_pop_vectors(F,trial_length,trial_idx,varargin)
% F is a long population vector from concatinated trials (of equal length trial_length)
% trial_idx saves indices of different trials types in F  
state_idx = [];
for v = 1:numel(varargin)
    if strcmpi(varargin{v},'state_idx')
        state_idx = varargin{v+1};
    end
end
pop_struct = struct();
fds = fields(trial_idx);
num_states = size(F,2);
if isempty(state_idx) % get certain states, otherwise take all
    state_idx = 1:num_states;
end

for f = 1:numel(fds)
    trial_num = trial_idx.(fds{f});
    this_F = nan(length(trial_num),trial_length,numel(state_idx));
    for i = 1:length(trial_num)
        this_F(i,:,:) = F((trial_num(i)-1)*(trial_length)+[1:trial_length],state_idx);
    end
    this_F_avg = squeeze(mean(this_F,1))';
    pop_struct.(fds{f}) = this_F;
%     pop_struct.([fds{f} '_avg']) = this_F_avg;
end

end

