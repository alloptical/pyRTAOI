% Adam Packer 20150319

% User set parameters
Filename = 'C:/Lloyd/20150319/ETLStairCase_0Vto5V_5Steps_StepEvery32-6ms_100reps.dat'; % Filename of output file, extension .dat
BiDi = 1; % bidirectional staircase
MinVolt = 0; % min voltage to output
MaxVolt = 5; % max voltage to output
NumStairs = 5; % total number of stairs including top and bottom stairs
TimeOnStair = 0.016300*2; % seconds
SampleRate = 10000; % samples per second
StaircaseReps = 1000; % repetitions of entire staircase pattern

if BiDi % deal with nonbidi some other day...
    % Calculated parameters
    NumStairsBiDi = (2 * NumStairs) - 2;
    SamplesInOneStair = round(TimeOnStair * SampleRate);
    SamplesInOneStairCase = round(NumStairsBiDi * TimeOnStair * SampleRate);
    StepsBetweenStairs = NumStairs - 1;
    VoltIncrementBetweenStairs = (MaxVolt - MinVolt) / StepsBetweenStairs;
    
    % Build a single staircase
    SingleStaircaseVector = zeros(1,SamplesInOneStairCase);
    AllStepValues = [MinVolt:VoltIncrementBetweenStairs:MaxVolt...
        fliplr(MinVolt+VoltIncrementBetweenStairs:VoltIncrementBetweenStairs:MaxVolt-VoltIncrementBetweenStairs)];
    for i = 1:NumStairsBiDi
        StartIdx = (i-1)*SamplesInOneStair+1;
        StopIdx = i*SamplesInOneStair;
        DesiredValue = AllStepValues(i);
        SingleStaircaseVector(StartIdx:StopIdx) = DesiredValue;
    end
    
    % Replicate staircases for final output
    Output = repmat(SingleStaircaseVector,1,StaircaseReps);
    Output(end+1) = MinVolt;
    
    % Check durations
    SingleRepDurationSecs = length(SingleStaircaseVector)/SampleRate;
    TotalDurationSecs = length(Output)/SampleRate;
end

% Visualise
figure
subplot(2,1,1)
plot((1:length(SingleStaircaseVector))/SampleRate,SingleStaircaseVector)
title('Single repetition')
xlabel('Seconds')
ylabel('Volts')
subplot(2,1,2)
plot((1:length(Output))/SampleRate,Output)
title('All repetitions')
xlabel('Seconds')
ylabel('Volts')

% Save
fid = fopen([Filename],'w','l');
fwrite(fid,Output,'double');
fclose(fid);
