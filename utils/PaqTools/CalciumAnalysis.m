function ExpDat=CalciumAnalysis(ExpDat,ThingsToAnalyze)

% Load movie timing file (txt, tab delimited, with one header row)
% if it is not there already
if ~isfield(ExpDat,'FluorescenceTime')
    MovieSyncs=pt_continuousabove(ExpDat.chanDev1_ai2_CameraSync,0,3,10,100000,0);
    if MovieSyncs(1,1)==1
        MovieStartIdx=MovieSyncs(6,1);
    else
        MovieStartIdx=MovieSyncs(5,1);
    end
    MovieStartTime=MovieStartIdx/10000;
    [filename,pathname,FilterIndex]=uigetfile('*.txt','Choose a TXT file with movie timing');
    if FilterIndex
        timingfile=fullfile(pathname,filename);
    else
        return
    end
    movietiming=dlmread(timingfile,'\t',1,0);
    ExpDat.FluorescenceTime=movietiming(:,3:4) + MovieStartTime;
end

figure;
Baseline=zeros(1,ThingsToAnalyze);
Noise=zeros(1,ThingsToAnalyze);
Signal=zeros(1,ThingsToAnalyze);
PercChange=zeros(1,ThingsToAnalyze);
SNR=zeros(1,ThingsToAnalyze);

for i=1:ThingsToAnalyze;
    subplot(ThingsToAnalyze,2,2*i - 1);
    EPidx=ExpDat.Analysis.FirstSpikesIdx(i);
    plot(ExpDat.chanDev1_ai0_VoltageCh1(EPidx-500:EPidx+1999));
    subplot(ThingsToAnalyze,2,2*i);
    idx=find(ExpDat.FluorescenceTime>ExpDat.Analysis.FirstSpikesIdx(i)/10000,1,'first');
    plot(ExpDat.Fluorescence(idx-25:idx+50));
    Baseline(i)=mean(ExpDat.Fluorescence(idx-25:idx));
    Noise(i)=std(ExpDat.Fluorescence(idx-25:idx));
    Signal(i)=max(ExpDat.Fluorescence(idx:idx+50)) - Baseline(i);
    PercChange(i)=(Signal(i)/Baseline(i)) * 100;
    SNR(i)=Signal(i)/Noise(i);
end;


ExpDat.Analysis.Baseline=Baseline;
ExpDat.Analysis.Noise=Noise;
ExpDat.Analysis.PercChange=PercChange;
ExpDat.Analysis.SNR=SNR;
ExpDat.Analysis.Signal=Signal;