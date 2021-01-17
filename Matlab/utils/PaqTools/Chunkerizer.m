%% Set us up the bomb
function varargout=Chunkerizer(varargin)

% Input
% 1) Incoming data to chunkerize
% Output
% 1) Chunks
% 2) Pulse Indices
% 3) Params

%% Bring it on
data=varargin{1};
try
    Params=varargin{2};
    if ~iscell(Params)
        Params={'1','2','0','0.1','10','10000','3','100','1','10','2000','0'};
    end
catch
    Params={'1','2','0','0.1','10','10000','3','100','1','10','500','0'};
end    
answer=inputdlg({'Recording Channel' 'Stimulation Channel' 'Pulse Baseline' 'Pulse Threshold' 'Pulse minimum time' ...
    'Pulse maximum time' 'Pulses per target' 'Targets (ignored if iterations=1)' ...
    'Number of iterations' 'Samples to get before each chunk', 'Samples to get after each chunk', ...
    'Milliseconds After Trigger To Ignore New Triggers'}, ...
    'Specify parameters',1,Params);
if isempty(answer)
    varargout{1}=Inf;
    varargout{2}=Inf;
    varargout{3}=Inf;
    varargout{4}=Inf;
    return
end
Params=struct('RecordChannel',str2double(answer(1)),'StimChannel',str2double(answer(2)),'Baseline',str2double(answer(3)), ...
    'Threshold',str2double(answer(4)),'Mintime',str2double(answer(5)),'Maxtime',str2double(answer(6)), ...
    'PulsesPerTarget',str2double(answer(7)),'NumTargets',str2double(answer(8)), ...
    'NumIterations',str2double(answer(9)),'PreChunkSamples',str2double(answer(10)), ...
    'PostChunkSamples',str2double(answer(11)),'Ignorepts',str2double(answer(12)));

%% Find pulse indices and chunkerize!
PulseIndices=pt_continuousabove(data(:,Params.StimChannel),Params.Baseline,Params.Threshold,Params.Mintime,Params.Maxtime,Params.Ignorepts);
if Params.NumIterations==1
    % Handle instances in which data is missing
    if mod(size(PulseIndices,1),Params.PulsesPerTarget)~=0
        q1=sprintf('%d pulses found in the data. \n',size(PulseIndices,1));
        q2=sprintf('Supposedly, there are %d pulses per target. \n',Params.PulsesPerTarget);
        q3=sprintf('Click yes to assume data is missing. \n');
        q4=sprintf('Click no to enter a new pulses per target value. \n');
        questans=questdlg([q1 q2 q3 q4],'Pulses found not divisible by pulses per target');
        if strcmp(questans,'Yes')
            DataMissing=1;
            Params.NumTargets=floor(size(PulseIndices,1)/Params.PulsesPerTarget);
        elseif strcmp(questans,'No')
            newanswer=APinputdlg({'Pulses Per Target'},'New value',1,{'1'});
            Params.PulsesPerTarget=str2double(newanswer(1));
            DataMissing=0;
        else
            varargout{1}=Inf;
            varargout{2}=Inf;
            varargout{3}=Inf;
            varargout{4}=Inf;
            return
        end
    else
        Params.NumTargets=size(PulseIndices,1)/Params.PulsesPerTarget;
        DataMissing=0;
    end
else
    DataMissing=NaN; % haven't written the code for more than 1 iteration
end

ChunkerizedDataIndices=zeros(Params.NumTargets,2,Params.NumIterations);
for i=0:Params.NumIterations-1
    for j=0:Params.NumTargets-1
        ThisChunkStartIdx=i*Params.NumTargets + j + 1;
        ThisChunkFirstPulseIdx=1+(Params.PulsesPerTarget*(ThisChunkStartIdx-1));
        ThisChunkStartSample=PulseIndices(ThisChunkFirstPulseIdx,1)-Params.PreChunkSamples;
        ThisChunkStopSample=PulseIndices(ThisChunkFirstPulseIdx+Params.PulsesPerTarget-1,2)+Params.PostChunkSamples;
        ChunkerizedDataIndices(j+1,1,i+1)=ThisChunkStartSample;
        ChunkerizedDataIndices(j+1,2,i+1)=ThisChunkStopSample;
    end
end
        
%% Send output
varargout{1}=ChunkerizedDataIndices;
varargout{2}=PulseIndices;
varargout{3}=Params;
varargout{4}=DataMissing;