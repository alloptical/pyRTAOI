function [output] = plot_psycho_curve(trials,plot_stim_types,plot_var_types)
%% set up stim matrix
% this should be identical to og_CrossStim_matrix_sinusoid
% need to enter these value manually from piezo calibration
R_1 = 0; % in volts
R_2 = 0.625;
R_3 = 1.25;
R_4 = 1.875;
R_5 = 2.5;

L_1 = 0;
L_2 = 0.625;
L_3 = 1.25;
L_4 = 1.875;
L_5 = 2.5;

LeftWhisker = [L_1,L_1,L_1,L_1,L_1;
               L_2,L_2,L_2,L_2,L_2;
               L_3,L_3,L_3,L_3,L_3;
               L_4,L_4,L_4,L_4,L_4;
               L_5,L_5,L_5,L_5,L_5];
           
RightWhisker = [R_1,R_2,R_3,R_4,R_5;
                R_1,R_2,R_3,R_4,R_5;
                R_1,R_2,R_3,R_4,R_5;
                R_1,R_2,R_3,R_4,R_5;
                R_1,R_2,R_3,R_4,R_5];
% vectorise matrix
RightWhisker = RightWhisker(:)';
LeftWhisker = LeftWhisker(:)';

idx_match.stim1 = [2,3,4,5,8,9,10,14,15,20];
idx_match.stim2 = [6,11,12,16,17,18,21,22,23,24];
idx_match.stim3 = [1,7,13,19,25];

%% get psychometric curve
num_vars = length(plot_var_types);
rightchoices = nan(1,num_vars);
leftchoices = nan(1,num_vars);
leftvolts =  nan(1,num_vars);
rightvolts =  nan(1,num_vars);
pc_rightchoices = nan(1,num_vars);
pc_leftchoices = nan(1,num_vars);
pc_misses = nan(1,num_vars);
for v = 1:num_vars
    this_var = plot_var_types(v);
    this_stim = plot_stim_types(v);    
    this_num_trials = numel(find(trials.stim_var == this_var & trials.stim_type == this_stim));
    this_idx = idx_match.(['stim' num2str(plot_stim_types(v))])(plot_var_types(v)+1);
    
    leftvolts(v) = LeftWhisker(this_idx);
    rightvolts(v) = RightWhisker(this_idx);
    
    rightchoices(v) = numel( find(trials.stim_var == this_var & trials.stim_type == this_stim & trials.firstresponse == 2));
    leftchoices(v) = numel( find(trials.stim_var == this_var & trials.stim_type == this_stim & trials.firstresponse == 1));   
    
    pc_rightchoices(v) = rightchoices(v)/this_num_trials;
    pc_leftchoices(v) = leftchoices(v)/this_num_trials;
    pc_misses(v) = numel(find(trials.stim_var == this_var & trials.stim_type == this_stim & trials.miss == 1))/this_num_trials;
end

sanity_check = double(pc_rightchoices+pc_leftchoices+pc_misses)

output = struct();
output.leftvolts = leftvolts;
output.rightvolts = rightvolts;
output.rightchoices = rightchoices;
output.pc_rightchoices = pc_rightchoices;
output.pc_misses = pc_misses;
output.pc_leftchoices = pc_leftchoices;
outpot.leftchoices = leftchoices;



end

