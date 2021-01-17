function [leftvolt,rightvolt] = get_volts(this_stim,this_var)
% return deflection voltages given stim type and var type
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
leftvolt = 0;
rightvolt = 0;
if this_stim<=3
    this_idx = idx_match.(['stim' num2str(this_stim)])(this_var+1);
    
    leftvolt = LeftWhisker(this_idx);
    rightvolt= RightWhisker(this_idx);
end
if this_var ==9 % easy trials with long deflection durations
    if this_stim==1
        leftvolt= L_5; rightvolt =0;
    elseif  this_stim==2
        leftvolt= 0; rightvolt =R_5;
    end
end

end

