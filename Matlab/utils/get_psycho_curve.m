function [output,num_vars] = get_psycho_curve(trials,stim_types,stim_vars)
%% get psychometric curve
num_vars = length(stim_vars);
rightchoices = nan(1,num_vars);
leftchoices = nan(1,num_vars);
leftvolts =  zeros(1,num_vars);
rightvolts =  zeros(1,num_vars);
pc_rightchoices = nan(1,num_vars);
pc_leftchoices = nan(1,num_vars);
pc_misses = nan(1,num_vars);
for v = 1:num_vars
    this_var = stim_vars(v);
    this_stim = stim_types(v);    
    this_num_trials = numel(find(trials.stim_var == this_var & trials.stim_type == this_stim));
    
    
%     if this_stim<=3 
%         this_idx = idx_match.(['stim' num2str(this_stim)])(this_var+1);
%         
%         leftvolts(v) = LeftWhisker(this_idx);
%         rightvolts(v) = RightWhisker(this_idx);
%     end
%     if this_var ==9 % easy trials with long deflection durations
%         if this_stim==1
%             leftvolts(v)= L_5; rightvolts(v) =0;
%         elseif  this_stim==2
%             leftvolts(v)= 0; rightvolts(v) =R_5;
%         end
%     end
   [leftvolts(v),rightvolts(v)] = get_volts(this_stim,this_var);
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
output.leftchoices = leftchoices;
output.stim_vars = stim_vars;
output.stim_types = stim_types;




end

