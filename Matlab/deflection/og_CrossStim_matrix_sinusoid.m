clear all

% Trial Sequence is:
% \Dropbox\Bruker1\Oliver\generate_trial_sequence\bruker1\og_CrossStim_trialseq.m

sample_rate_hz          = 1000;
num_diff_output_lines   = 4;

daq_session = daq.createSession('ni');
daq_session2 = daq.createSession('ni');
data = zeros(5000, num_diff_output_lines);

% stim channels
addDigitalChannel(daq_session, 'Dev1', 'port1/line2', 'InputOnly');
addDigitalChannel(daq_session, 'Dev1', 'port1/line3', 'InputOnly');
addDigitalChannel(daq_session, 'Dev1', 'port0/line24', 'InputOnly');
addDigitalChannel(daq_session, 'Dev1', 'port0/line25', 'InputOnly');
addDigitalChannel(daq_session, 'Dev1', 'port1/line0', 'InputOnly');
% variation channels
addDigitalChannel(daq_session, 'Dev1', 'port0/line11', 'InputOnly');
addDigitalChannel(daq_session, 'Dev1', 'port0/line12', 'InputOnly');
addDigitalChannel(daq_session, 'Dev1', 'port0/line13', 'InputOnly');
addDigitalChannel(daq_session, 'Dev1', 'port0/line14', 'InputOnly');

% analog output channels
addAnalogOutputChannel(daq_session2, 'Dev1','ao2', 'Voltage');
addAnalogOutputChannel(daq_session2, 'Dev1','ao1', 'Voltage');
addAnalogOutputChannel(daq_session2, 'Dev1','ao0', 'Voltage'); % for lick port motor
addAnalogOutputChannel(daq_session2, 'Dev1','ao3', 'Voltage'); % Photostim trig

lickportmotor = zeros(1,5000); lickportmotor(1330:3330) = 5;
% set up photostim trigger

photostim_delay = 1;% i.e. factor in sensory delay to cortex?
photostim_trig_v = 5; %Volt
photostim_duration = 20; %in samples, 

photostim_trig = zeros(5000,1);
photostim_trig(photostim_delay:photostim_delay+photostim_duration) = photostim_trig_v; 

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

% set up stim matrix
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

stim1 = [2,3,4,5,8,9,10,14,15,20];
stim2 = [6,11,12,16,17,18,21,22,23,24];
stim3 = [1,7,13,19,25];


% analog output channels
LeftPiezo = 1;
RightPiezo = 2;
LickportMotor = 4;
PhotostimTrig = 3;

sineFreq = 20; % hz
sineDur = 0.5; % seconds
sineAmp = 1; % set it to 1V here then scale it accordingli in the loop below

SineWave_short = sin(2*pi*sineFreq*[0:1/sample_rate_hz:sineDur-(1/sample_rate_hz)]) * sineAmp; 
SineWave_long = sin(2*pi*sineFreq*[0:1/sample_rate_hz:1.5-(1/sample_rate_hz)]) * sineAmp; 


%% baseline behaviour imaging
disp('Waiting for first trigger...')


while true
     
    scanTriggers = [daq_session.inputSingleScan];
    trigger_states = scanTriggers(1:5);
    VarStates = scanTriggers(6:9);   
    
    if any(trigger_states)                    
        which_stim = find(trigger_states);    
        varNum = bi2de(VarStates)+1;    
        if which_stim == 1
                data(1:500,LeftPiezo) = SineWave_short*LeftWhisker(stim1(varNum))';
                data(1:500,RightPiezo) = SineWave_short*RightWhisker(stim1(varNum))';
                if varNum == 10
                    data(1:1500,LeftPiezo) = SineWave_long*L_5';
                    data(1:1500,RightPiezo) = SineWave_long*0';
                end
        elseif which_stim == 2 % right > left stim
                data(1:500,LeftPiezo) = SineWave_short*LeftWhisker(stim2(varNum))';
                data(1:500,RightPiezo) = SineWave_short*RightWhisker(stim2(varNum))';
                if varNum == 10
                    data(1:1500,LeftPiezo) = SineWave_long*0';
                    data(1:1500,RightPiezo) = SineWave_long*R_5';
                end
        elseif which_stim == 3 % no stim catch
               data(1:500,LeftPiezo) = SineWave_short*LeftWhisker(stim3(varNum))';
               data(1:500,RightPiezo) = SineWave_short*RightWhisker(stim3(varNum))';
        end
               
        data(:,LickportMotor) = lickportmotor; % motor is driven on every trial  
        
        queueOutputData(daq_session2, data);
        daq_session2.startForeground();        
              
        disp(['Characterisation Imaging. Trial delivered: Stim ' num2str(which_stim) ', Variation ' num2str(varNum)])
        data = zeros(5000, num_diff_output_lines);
    end
end


%% STIM EXPERIMENT

% manually change these values
LeftPiezoThresh = 1.25;
RightPiezoThresh = 1.25;

while true
     
    scanTriggers = [daq_session.inputSingleScan];
    trigger_states = scanTriggers(1:5);
    VarStates = scanTriggers(6:9);   
    
    if any(trigger_states)                    
        which_stim = find(trigger_states);    
        varNum = bi2de(VarStates)+1;    

        if which_stim == 1 % left > right stim
            if varNum == 1
                data(1:500,LeftPiezo) = SineWave_short*L_5';
                data(1:500,RightPiezo) = SineWave_short*0';
            elseif varNum == 2
                data(1:500,LeftPiezo) = SineWave_short*LeftPiezoThresh';   % PHOTOSTIM + THRESHOLD STIMULUS
                data(1:500,RightPiezo) = SineWave_short*RightPiezoThresh';                
            elseif varNum == 3
                data(1:1500,LeftPiezo) = SineWave_long*L_5';   % PHOTOSTIM + THRESHOLD STIMULUS
                data(1:1500,RightPiezo) = SineWave_long*0';                
            end
        elseif which_stim == 2 % right > left stim
            if varNum == 1
                data(1:500,LeftPiezo) = SineWave_short*0';
                data(1:500,RightPiezo) = SineWave_short*R_5';
            elseif varNum == 2
                data(1:500,LeftPiezo) = SineWave_short*LeftPiezoThresh';   % PHOTOSTIM + THRESHOLD STIMULUS
                data(1:500,RightPiezo) = SineWave_short*RightPiezoThresh';                
           elseif varNum == 3
                data(1:1500,LeftPiezo) = SineWave_long*0';   % PHOTOSTIM + THRESHOLD STIMULUS
                data(1:1500,RightPiezo) = SineWave_long*R_5';                
            end
        elseif which_stim == 3 % no stim catch
               data(1:500,LeftPiezo) = SineWave_short*LeftWhisker(stim3(varNum))';
                data(1:500,RightPiezo) = SineWave_short*RightWhisker(stim3(varNum))';
                
        elseif which_stim == 4 % sensory + photostim trials 
            if varNum == 0
            	data(1:500,LeftPiezo) = 0;   % PHOTOTSTIM ALONE, NO SENSORY
                data(1:500,RightPiezo) =0;                
            else 
                data(1:500,LeftPiezo) = SineWave_short*LeftPiezoThresh';   % PHOTOSTIM + THRESHOLD STIMULUS
                data(1:500,RightPiezo) = SineWave_short*RightPiezoThresh';                
            end
           data(:,PhotostimTrig) = photostim_trig;  
                     
        elseif which_stim == 5 % photostim trials only
       
           data(:,PhotostimTrig) = photostim_trig;  
        end
               
        data(:,LickportMotor) = lickportmotor; % motor is driven on every trial  
        
        queueOutputData(daq_session2, data);
        daq_session2.startForeground();        
              
        disp(['Stimulation Experiment. Trial delivered: Stim ' num2str(which_stim) ', Variation ' num2str(varNum)])

        data = zeros(5000, num_diff_output_lines);
    end
end

%% widefield 

flipflop = repmat([0.5,1.5,2.5],1,50)
count = 1;
while true
     if count > 4
         count = 1;
     end
    scanTriggers = [daq_session.inputSingleScan];
    trigger_states = scanTriggers(1:5)
    VarStates = scanTriggers(6:9);   
    
    if any(trigger_states)                    
        which_stim = find(trigger_states)       
        varNum = bi2de(VarStates)+1;        
        if which_stim == 1 % right > left stim
                data(1:500,LeftPiezo) = SineWave_short*4';
                data(1:500,RightPiezo) = SineWave_short*0';
                count = count + 1;
                
        elseif which_stim == 2 % right > left stim
               data(1:500,LeftPiezo) = SineWave_short*L_5';
                data(1:500,RightPiezo) = SineWave_short*0';
                data(100:end,LeftPiezo) = 0';

        end
                       
        queueOutputData(daq_session2, data);
        daq_session2.startForeground();        
              
        disp(['Trial delivered: Stim ' num2str(which_stim) ', Variation ' num2str(varNum)])

        data = zeros(5000, num_diff_output_lines);
    end
end


%% LED bias experiment
disp('Waiting for first trigger...')

leftAMP = [0 L_1 L_2 L_3 L_4 L_5]
rightAMP = [0 L_5 L_4 L_3 L_2 L_1]

while true
     
    scanTriggers = [daq_session.inputSingleScan];
    trigger_states = scanTriggers(1:5);
    VarStates = scanTriggers(6:9);   
    
    if any(trigger_states)                    
        which_stim = find(trigger_states);    
        varNum = bi2de(VarStates)+1;    
        if which_stim == 1
                data(1:500,LeftPiezo) = SineWave_short*LeftWhisker(stim1(varNum))';
                data(1:500,RightPiezo) = SineWave_short*RightWhisker(stim1(varNum))';
                if varNum == 10
                    data(1:1500,LeftPiezo) = SineWave_long*L_5';
                    data(1:1500,RightPiezo) = SineWave_long*0';
                end
        elseif which_stim == 2 % right > left stim
                data(1:500,LeftPiezo) = SineWave_short*LeftWhisker(stim2(varNum))';
                data(1:500,RightPiezo) = SineWave_short*RightWhisker(stim2(varNum))';
                if varNum == 10
                    data(1:1500,LeftPiezo) = SineWave_long*0';
                    data(1:1500,RightPiezo) = SineWave_long*R_5';
                end
        elseif which_stim == 3 % no stim catch
              data(1:500,LeftPiezo) = SineWave_short*LeftWhisker(stim3(varNum))';
               data(1:500,RightPiezo) = SineWave_short*RightWhisker(stim3(varNum))';
       elseif which_stim == 4 % no stim catch
               data(1:500,LeftPiezo) = SineWave_short*leftAMP((varNum))';
               data(1:500,RightPiezo) = SineWave_short*rightAMP((varNum))';
               data(1:500,3) = SineWave_short*5';

        end
               
        data(:,LickportMotor) = lickportmotor; % motor is driven on every trial  
        
        queueOutputData(daq_session2, data);
        daq_session2.startForeground();        
              
        disp(['Characterisation Imaging. Trial delivered: Stim ' num2str(which_stim) ', Variation ' num2str(varNum)])
        data = zeros(5000, num_diff_output_lines);
    end
end

