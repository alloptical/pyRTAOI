% adapted from experiment_trial_seq from Oli
clear variables

x = inputdlg({'ID'}, 'Enter Mouse ID',[1 50],{'OG'});

num_loops = 30;

%% For baseline session

% stim 1 - left discrimination trials (10 conditions)
% stim 2 - right discrimination trials (10 conditions)
% stim 3 - threshold trials (+ nogo catch)
% stim 4 - sensory + photostim (10 conditions)
% stim 5 - motor only (no cue)
% stim 6 - cue only (no motor)

stims = cell(1,num_loops);
vars = cell(1,num_loops);
photostim_var = cell(1,num_loops);

for loop = 1:num_loops

    if loop == 1
        stims{loop} = [1,2,1,2,1,2,1,2,1,2];%2,1,2,1,2,1,2,1,2]; % first set of trials are 100% easy trials in interleaved order
        vars{loop} =  [9 9 9 9 9 9 9 9 9 9];%,6,3,6,3,6,3,6,3,6];      
    else
        % baseline imaging
       stims{loop} =  [1 1 2 2 3 1 1 2 2 3 5 6]; % 20190706 - changed this slightly so its just the diagonal psychometric curve
       vars{loop} =   [3 5 4 6 2 3 5 4 6 2 0 0];   
       
   % passive square
     % stims{loop} = [1 1 1 2 2 2 3 3 3];
      %   vars{loop} =  [1,3,8,1,8,6,0,2,4];   
        photostim_var_rand = randperm(10);
        photostim_var{loop} = 0:9;
        photostim_var{loop} = photostim_var{loop}(photostim_var_rand);
        
        while true
            rand_ind = randperm(length(stims{loop}));
            stimtemp = stims{loop}(rand_ind);
            if all(abs((diff(diff(stimtemp))))>0) % diff diff checks for > 3 in a row
                vars{loop} = vars{loop}(rand_ind);
                stims{loop} = stimtemp;
                break
            end                
        end    
    end
end

trial_seq(1,:) = horzcat(stims{:});
trial_seq(2,:) = horzcat(vars{:});



easyLeftTrial = 10:6:10000;
easyRightTrial = 13:6:10000;

stim_count = 1;

try
    for i = 1:900
        if any(easyLeftTrial== i) 
            final_stim(i) = 1; 
            final_var(i) = 9;%easyLeftvar(c1);c1 = c1+1;
        elseif any(easyRightTrial== i) 
            final_stim(i) = 2; 
            final_var(i) = 9;%easyRightvar(c2);c2=c2+1;
        else
            final_stim(i) = trial_seq(1,stim_count); 
            final_var(i) = trial_seq(2,stim_count); stim_count = stim_count +1;
        end
    end 
catch
end

clear trial_seq

trial_seq(1,:) = final_stim;
trial_seq(2,:) = final_var;

%
fname = [x{1} '_baseline_imaging_experiment_' date ];
disp('Select directory to save trial sequence')
pause(1)
save_dir = uigetdir();
dlmwrite([save_dir filesep fname '.txt'],trial_seq)



%% PHOTOSTIMULATION EXPERIMENT TRIAL SEQUENCE
clear trial_seq stim var

clear variables

x = inputdlg({'ID'}, 'Enter Mouse ID',[1 50],{'OG'});
num_loops = 20;

num_loops = 70;

stims = cell(1,num_loops);
vars = cell(1,num_loops);
photostim_var = cell(1,num_loops);

for loop = 1:num_loops

    if loop == 1
        stims{loop} = [1,2,1,2,1,2,1,2,1,2];%2,1,2,1,2,1,2,1,2]; % first set of trials are 100% easy trials in interleaved order
        vars{loop} =  [2 2 2 2 2 2 2 2 2 2];%,6,3,6,3,6,3,6,3,6];      
    else              
     stims{loop} =  [1 1 2 2 3 4 4 5 5];
       vars{loop} = [1 0 0 1 0 0 1 0 1];   
       
        photostim_var_rand = randperm(10);
        photostim_var{loop} = 0:9;
        photostim_var{loop} = photostim_var{loop}(photostim_var_rand);
        
%         while true
            rand_ind = randperm(length(stims{loop}));
            stimtemp = stims{loop}(rand_ind);
%             if all(abs((diff(diff(stimtemp))))>0) % diff diff checks for > 3 in a row
                vars{loop} = vars{loop}(rand_ind);
                stims{loop} = stimtemp;
%                 break
%             end                
%         end    
    end
end

trial_seq(1,:) = horzcat(stims{:});
trial_seq(2,:) = horzcat(vars{:});



easyLeftTrial = 10:6:10000;
easyRightTrial = 13:6:10000;

stim_count = 1;

c1 = 1; c2 = 1;
try
    for i = 1:900
        if any(easyLeftTrial== i) 
            final_stim(i) = 1; 
            final_var(i) = 2;%easyLeftvar(c1);c1 = c1+1;
        elseif any(easyRightTrial== i) 
            final_stim(i) = 2; 
            final_var(i) = 2;%easyRightvar(c2);c2=c2+1;
        else
            final_stim(i) = trial_seq(1,stim_count); 
            final_var(i) = trial_seq(2,stim_count); stim_count = stim_count +1;
        end
    end 
catch
end

clear trial_seq


trial_seq(1,:) = final_stim;
trial_seq(2,:) = final_var;


fname = [x{1} '_stim_experiment_' date ];
disp('Select directory to save trial sequence')
pause(1)
save_dir = uigetdir();
dlmwrite([save_dir filesep fname '.txt'],trial_seq)


%

TrialOrder = trial_seq;
NaparmFolder = uigetdir();

ChannelsWithLaser = [4 5];

PMFolder = [NaparmFolder filesep 'PhaseMasks'];
BehavPMFolder = [NaparmFolder filesep 'BehavPhaseMasks'];
mkdir(BehavPMFolder);
PMs = glob([PMFolder filesep '*.tiff']);
xmlFiles = glob([NaparmFolder filesep '*.xml']);

%load the original naparm xml
fileID = fopen(xmlFiles{1},'r');
xml = fscanf(fileID, '%c');
fclose(fileID);
elementStarts = strfind(xml, '<PVMarkPointElement ');
elementStops = strfind(xml, '</PVMarkPointElement>') + numel('</PVMarkPointElement>') -1;
elements = [];
for i = 1:numel(elementStarts)
    elements{i} = xml(elementStarts(i) : elementStops(i));
end
behavXML = xml(1:elementStarts(2)-1);  % to have dummy included

% grow the xml file and copy phase mask to new folder in correct order
runningCount = 0;
for t = 1:size(TrialOrder,2)
    thisStim = TrialOrder(1,t);
    thisVar = TrialOrder(2,t);
    if any(thisStim==ChannelsWithLaser)
        runningCount = runningCount+1;
        fname = [num2str(runningCount, '%03i') '_Stim' num2str(thisStim) '_Var' num2str(thisVar) '_Trial' num2str(t, '%03i') '.tiff'];
        disp(fname)
        copyfile(PMs{thisVar+1}, [BehavPMFolder filesep fname])
        behavXML = [behavXML elements{thisVar+1+1}];  % plus 2 because 0-indexing and dummy
    end
end
behavXML = [behavXML xml(elementStops(end)+1:end)];

% save the new XML file
fileID = fopen(strrep(xmlFiles{1}, '.xml', 'Behav.xml'), 'w');
fprintf(fileID, '%s', behavXML);
fclose(fileID);

disp('Done')


%% PASSIVE SQUARE

clear variables

x = inputdlg({'ID'}, 'Enter Mouse ID',[1 50],{'OG'});

num_loops = 30;

%

% stim 1 - left discrimination trials (10 conditions)
% stim 2 - right discrimination trials (10 conditions)
% stim 3 - threshold trials (+ nogo catch)
% stim 4 - sensory + photostim (10 conditions)
% stim 5 - motor only (no cue)
% stim 6 - cue only (no motor)

stims = cell(1,num_loops);
vars = cell(1,num_loops);
photostim_var = cell(1,num_loops);

for loop = 1:num_loops

    if loop == 100
        stims{loop} = [1,2,1,2,1,2,1,2,1,2];%2,1,2,1,2,1,2,1,2]; % first set of trials are 100% easy trials in interleaved order
        vars{loop} =  [9 9 9 9 9 9 9 9 9 9];%,6,3,6,3,6,3,6,3,6];      
    else
        % baseline imaging

       
   % passive square
     stims{loop} = [1 1 1 2 2 2 3 3 3];
         vars{loop} =  [1,3,8,1,8,6,0,2,4];   
        photostim_var_rand = randperm(10);
        photostim_var{loop} = 0:9;
        photostim_var{loop} = photostim_var{loop}(photostim_var_rand);
        
        while true
            rand_ind = randperm(length(stims{loop}));
            stimtemp = stims{loop}(rand_ind);
            if all(abs((diff(diff(stimtemp))))>0) % diff diff checks for > 3 in a row
                vars{loop} = vars{loop}(rand_ind);
                stims{loop} = stimtemp;
                break
            end                
        end    
    end
end

trial_seq(1,:) = horzcat(stims{:});
trial_seq(2,:) = horzcat(vars{:});



easyLeftTrial = 300:6:10000;
easyRightTrial = 303:6:10000;

stim_count = 1;

try
    for i = 1:900
        if any(easyLeftTrial== i) 
            final_stim(i) = 1; 
            final_var(i) = 9;%easyLeftvar(c1);c1 = c1+1;
        elseif any(easyRightTrial== i) 
            final_stim(i) = 2; 
            final_var(i) = 9;%easyRightvar(c2);c2=c2+1;
        else
            final_stim(i) = trial_seq(1,stim_count); 
            final_var(i) = trial_seq(2,stim_count); stim_count = stim_count +1;
        end
    end 
catch
end

clear trial_seq

trial_seq(1,:) = final_stim;
trial_seq(2,:) = final_var;

%
fname = [x{1} '_passive_square_imaging_experiment_' date ];
disp('Select directory to save trial sequence')
pause(1)
save_dir = uigetdir();
dlmwrite([save_dir filesep fname '.txt'],trial_seq)


%%



clear variables

x = inputdlg({'ID'}, 'Enter Mouse ID',[1 50],{'OG'});

num_loops = 30;

%

% stim 1 - left discrimination trials (10 conditions)
% stim 2 - right discrimination trials (10 conditions)
% stim 3 - threshold trials (+ nogo catch)
% stim 4 - sensory + photostim (10 conditions)
% stim 5 - motor only (no cue)
% stim 6 - cue only (no motor)

stims = cell(1,num_loops);
vars = cell(1,num_loops);
photostim_var = cell(1,num_loops);

for loop = 1:num_loops

    if loop == 1
        stims{loop} = [1,2,1,2,1,2,1,2,1,2];%2,1,2,1,2,1,2,1,2]; % first set of trials are 100% easy trials in interleaved order
        vars{loop} =  [9 9 9 9 9 9 9 9 9 9];%,6,3,6,3,6,3,6,3,6];      
    else
        % baseline imaging
       stims{loop} =  [1 1 1 2 2 2 3 3 1 1 2 2 3 4 4 4 4 4 4]; % 20190706 - changed this slightly so its just the diagonal psychometric curve
       vars{loop} =   [3 3 5 4 6 6 2 2 3 5 4 6 0 0 1 2 3 4 5];   
       
   % passive square
     % stims{loop} = [1 1 1 2 2 2 3 3 3];
      %   vars{loop} =  [1,3,8,1,8,6,0,2,4];   
        photostim_var_rand = randperm(10);
        photostim_var{loop} = 0:9;
        photostim_var{loop} = photostim_var{loop}(photostim_var_rand);
        
        while true
            rand_ind = randperm(length(stims{loop}));
            stimtemp = stims{loop}(rand_ind);
            if all(abs((diff(diff(stimtemp))))>0) % diff diff checks for > 3 in a row
                vars{loop} = vars{loop}(rand_ind);
                stims{loop} = stimtemp;
                break
            end                
        end    
    end
end

trial_seq(1,:) = horzcat(stims{:});
trial_seq(2,:) = horzcat(vars{:});



easyLeftTrial = 10:6:10000;
easyRightTrial = 13:6:10000;

stim_count = 1;

try
    for i = 1:900
        if any(easyLeftTrial== i) 
            final_stim(i) = 1; 
            final_var(i) = 9;%easyLeftvar(c1);c1 = c1+1;
        elseif any(easyRightTrial== i) 
            final_stim(i) = 2; 
            final_var(i) = 9;%easyRightvar(c2);c2=c2+1;
        else
            final_stim(i) = trial_seq(1,stim_count); 
            final_var(i) = trial_seq(2,stim_count); stim_count = stim_count +1;
        end
    end 
catch
end

clear trial_seq

trial_seq(1,:) = final_stim;
trial_seq(2,:) = final_var;

%
fname = [x{1} '_LEDbiasexperiment_' date ];
disp('Select directory to save trial sequence')
pause(1)
save_dir = uigetdir();
dlmwrite([save_dir filesep fname '.txt'],trial_seq)
