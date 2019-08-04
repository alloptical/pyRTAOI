function [output_shuff_trial_order,output_shuff_trial_type,final_max_consec] = get_shuff_trial_idx(tex_trial_types_idx,max_num_consecutive)
max_num_shuff = 20;
num_trials = length(tex_trial_types_idx);
this_max = num_trials;
num_shuff = 0;
all_trial_types = unique(tex_trial_types_idx);
num_trial_types = length(all_trial_types);

while num_shuff<max_num_shuff && this_max> max_num_consecutive
    this_shuff_trial_order = randperm(num_trials);
    this_shuff_trial_type = tex_trial_types_idx(this_shuff_trial_order);
    temp_max = nan(1,num_trial_types);
    for tt = 1:num_trial_types
        this_trial_type = all_trial_types(tt);
        temp = find(this_shuff_trial_type==this_trial_type);
        temp = [temp(1) find(this_shuff_trial_type==this_trial_type) temp(end)];
        
        % find longes consecutive trial type in array
        D = diff(diff(temp)==1);
        
        % make the first non-zero D positive
        temp_idx = find(D~=0);
        if D(temp_idx(1))<0
            D = -D;
        end
        %
        B = find([D>0]);
        E = find([D<0,true])+1;
        
        if length(E)>length(B)
            E = E(1:length(B));
        end
        
        if length(B)>length(E)
            B = B(2:length(B));
        end
        
        [temp_max(tt),idx] = max(E-B);        
        
%         find(this_shuff_trial_type==this_trial_type)
%         D
%         B
%         E
    end
    
    [this_max, this_type] = max(temp_max);
    num_shuff = num_shuff+1;
end

output_shuff_trial_order = this_shuff_trial_order;
output_shuff_trial_type= this_shuff_trial_type;
final_max_consec = this_max;

disp(['max number consecutive: ' num2str(final_max_consec) ' for type ' num2str(all_trial_types(this_type)) ' number shuffles: ' num2str(num_shuff)])

end

