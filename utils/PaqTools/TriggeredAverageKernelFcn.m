function outputfig = TriggeredAverageKernelFcn(trigtraces, rawtraces, zeroedtraces, SampleRate, trigtimes);
%varargin should contain
% - trigtraces
% - rawtraces
% - zeroedtraces

meantrigtraces = mean(trigtraces,1);
meanrawtraces = mean(rawtraces,1);
meanzeroedtraces = mean(zeroedtraces,1);
msecpersample= 1000 / SampleRate;

data.data{1} = trigtraces;
data.data{2} = meantrigtraces;
data.data{3} = rawtraces;
data.data{4} = meanrawtraces;
data.data{5} = zeroedtraces;
data.data{6} = meanzeroedtraces;
data.data{7} = trigtimes;
data.names{1} = 'Triggers';
data.names{2} = 'AveragedTriggers';
data.names{3} = 'TriggeredTraces';
data.names{4} = 'AveragedTraces';
data.names{5} = 'ZeroedTraces';
data.names{6} = 'AveragedZeroedTraces';
data.names{7} = 'TriggerTimes';

if isempty(trigtraces);
    errordlg('No Triggers Detected')
    return
end

outputfig = figure('units','normalized','position',[0 .05 .3 .85]);
TrigAvgExportMenu = uimenu('Label','Export','Parent',outputfig);
uimenu(TrigAvgExportMenu(end),'Label','To Base Workspace','callback',@ExportPeakTrigAvgToBaseWorkspace);
setappdata(outputfig,'data',data)

subplot(6,1,1);
lastx=(size(trigtraces,2)-1)*msecpersample;
x=0:msecpersample:lastx;
plot(x,trigtraces');
xlim([0 lastx])
maxval = max(max(trigtraces));
minval = min(min(trigtraces));
margin = .1*(maxval - minval);
ylim([minval-margin maxval+margin]);
title(['Triggers: ',num2str(size(trigtraces,1)),' total.'])

subplot(6,1,2);
lastx=(size(meantrigtraces,2)-1)*msecpersample;
x=0:msecpersample:lastx;
plot(x,meantrigtraces');
xlim([0 lastx])
maxval = max(max(meantrigtraces));
minval = min(min(meantrigtraces));
margin = .1*(maxval - minval);
ylim([minval-margin maxval+margin]);
title ('Averaged Triggers')

subplot(6,1,3);
lastx=(size(rawtraces,2)-1)*msecpersample;
x=0:msecpersample:lastx;
plot(x,rawtraces');
xlim([0 lastx])
maxval = max(max(rawtraces));
minval = min(min(rawtraces));
margin = .1*(maxval - minval);
ylim([minval-margin maxval+margin]);
title ('Triggered Traces')

subplot(6,1,4);
lastx=(size(meanrawtraces,2)-1)*msecpersample;
x=0:msecpersample:lastx;
plot(x,meanrawtraces');
xlim([0 lastx])
maxval = max(max(meanrawtraces));
minval = min(min(meanrawtraces));
margin = .1*(maxval - minval);
ylim([minval-margin maxval+margin]);
title ('Averaged Triggered Traces')

subplot(6,1,5);
lastx=(size(zeroedtraces,2)-1)*msecpersample;
x=0:msecpersample:lastx;
plot(x,zeroedtraces');
xlim([0 lastx])
maxval = max(max(zeroedtraces));
minval = min(min(zeroedtraces));
margin = .1*(maxval - minval);
ylim([minval-margin maxval+margin]);
title ('Zeroed Triggered Traces.  (Zero set by last point before trigger)')

subplot(6,1,6);
lastx=(size(meanzeroedtraces,2)-1)*msecpersample;
x=0:msecpersample:lastx;
plot(x,meanzeroedtraces');
xlim([0 lastx])
maxval = max(max(meanzeroedtraces));
minval = min(min(meanzeroedtraces));
margin = .1*(maxval - minval);
ylim([minval-margin maxval+margin]);
title ('Averaged Zeroed Traces');