function [cell_sta_struct] = add_pool_trials(cell_sta_struct)
    % balance correct and incorrect number of trials
    stim_size_1 = min(size(cell_sta_struct(1).stim_1_correct,2),size(cell_sta_struct(1).stim_1_incorrect,2));
    stim_rand_correct_1 = randperm( size(cell_sta_struct(1).stim_1_correct,2),stim_size_1);
    stim_rand_incorrect_1 = randperm( size(cell_sta_struct(1).stim_1_incorrect,2),stim_size_1);
    
    stim_size_2 = min(size(cell_sta_struct(1).stim_2_correct,2),size(cell_sta_struct(1).stim_2_incorrect,2));
    stim_rand_correct_2 = randperm( size(cell_sta_struct(1).stim_2_correct,2),stim_size_2);
    stim_rand_incorrect_2 = randperm( size(cell_sta_struct(1).stim_2_incorrect,2),stim_size_2);
    
    port_size_1 = min(size(cell_sta_struct(1).stim_1_correct,2),size(cell_sta_struct(1).stim_2_incorrect,2));
    port_rand_correct_1 = randperm( size(cell_sta_struct(1).stim_1_correct,2),port_size_1);
    port_rand_incorrect_1 = randperm( size(cell_sta_struct(1).stim_2_incorrect,2),port_size_1);
    
    port_size_2 = min(size(cell_sta_struct(1).stim_2_correct,2),size(cell_sta_struct(1).stim_1_incorrect,2));
    port_rand_correct_2 = randperm( size(cell_sta_struct(1).stim_2_correct,2),port_size_2);
    port_rand_incorrect_2 = randperm( size(cell_sta_struct(1).stim_1_incorrect,2),port_size_2);
    num_cells = size(cell_sta_struct,2);
    for c = 1:num_cells
        cell_sta_struct(c).stim_1 = [cell_sta_struct(c).stim_1_correct(:,stim_rand_correct_1) cell_sta_struct(c).stim_1_incorrect(:,stim_rand_incorrect_1)];
        cell_sta_struct(c).stim_2 = [cell_sta_struct(c).stim_2_correct(:,stim_rand_correct_2) cell_sta_struct(c).stim_2_incorrect(:,stim_rand_incorrect_2)];
        
        % randomly select a subset of correct trials to match numbers of incorrect trials
        cell_sta_struct(c).port_1 = [cell_sta_struct(c).stim_1_correct(:,port_rand_correct_1) cell_sta_struct(c).stim_2_incorrect(:,port_rand_incorrect_1)];
        cell_sta_struct(c).port_2 = [cell_sta_struct(c).stim_2_correct(:,port_rand_correct_2) cell_sta_struct(c).stim_1_incorrect(:,port_rand_incorrect_2)];
        
    end

end

