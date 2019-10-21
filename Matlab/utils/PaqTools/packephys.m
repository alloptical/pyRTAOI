% Adam Packer
% February 11th, 2007
% Based on Toledo Rodriguez, Markram etc genetics predicts ephys
% (all of markram's e1-e61 params are included)

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Choose and open file
function [Ephys]=packephys(fullpath)
fid=fopen(fullpath);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read in rate, number of channels, and channel names
rate=fread(fid,1,'float32','b');
numchans=fread(fid,1,'float32','b');
for i=1:numchans;
    number_of_characters=fread(fid,1,'float32','b');
    channelname{i}=[];
    for j=1:number_of_characters
        channelname{i}=[channelname{i}, strrep(fread(fid,1,'float32=>char', 'b'),' ','')];
    end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read in hardware channel ('HWchan') and units if available (*.paq only)
[pathstr, name, ext, versn] = fileparts(fullpath);
if strcmp(ext,'.paq')
    for k=1:numchans
        number_of_characters=fread(fid,1,'float32','b');
        HWchan{k}=[];
        for m=1:number_of_characters
            HWchan{k}=[HWchan{k}, strrep(fread(fid,1,'float32=>char','b'),' ','')];
        end
    end
    for n=1:numchans
        number_of_characters=fread(fid,1,'float32','b');
        units{n}=[];
        for q=1:number_of_characters
            units{n}=[units{n}, strrep(fread(fid,1,'float32=>char','b'),' ','')];
        end
    end
elseif strcmp(ext,'.bin')
    for k=1:numchans
        HWchan{k}='unknown';
        units{k}='unknown';
    end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get filesize and fposition (which should be after all header info)
dirinfo=dir(pathstr);
for n=1:size(dirinfo);
    if strcmp(dirinfo(n).name,[name ext]);
        filesize=dirinfo(n).bytes;
    end
end
fposition=ftell(fid);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read in data 
fseek(fid,fposition,'bof');
data=fread(fid,[numchans,filesize/numchans],'*float32','b');
data=data';
fclose(fid);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup chunking data
% Multiply by 2 to transform from output to input sampling rate (5 kHz to
% 10 kHz) Also always take 500 ms on the front and back of each chunk for
% pre and post injection calculations.
% AP_Drop=2*(2500:0.5:9999.5);
% AP_Waveform=2*(10000:0.5:15249.5);
% IV=[];
% for a=0:10
%     IVstart=15250+2*a*5000;
%     IVstop=20249+(2*a+1)*5000;
%     IV(a+1,:)=2*[IVstart:0.5:IVstop+0.5];
% end
% InputR_HypPre=2*[125250:0.5:133584.5];
% InputR_HypTest=2*[128585:0.5:134584.5];
% InputR_HypPost=2*[129585:0.5:135249.5];
% InputR_DepPre=2*[135250:0.5:143584.5];
% InputR_DepTest=2*[138585:0.5:144584.5];
% InputR_DepPost=2*[139585:0.5:145249.5];
% Delta=2*[145250:0.5:150259.5];
% Ramp=2*[150260:0.5:157759.5];
% Discharge=[];
% for b=0:10
%     DischargeStart=157760+b*6000;
%     DischargeStop=162759+b*6000+1000;
%     Discharge(b+1,:)=2*[DischargeStart:0.5:DischargeStop+0.5];
% end

% The following are WITHOUT pre and post 500 ms
% AP_Drop=2*[5000:0.5:7499.5];
% AP_Waveform=2*2*[12500:0.5:12749.5];
% IV=[];
% for a=0:10
%     IVstart=17750+2*a*5000;
%     IVstop=17749+(2*a+1)*5000;
%     IV(a+1,:)=2*[IVstart:0.5:IVstop+0.5];
% end
% InputR_HypPre=2*[127750:0.5:131084.5];
% InputR_HypTest=2*[131085:0.5:132084.5];
% InputR_HypPost=2*[132085:0.5:132749.5];
% InputR_DepPre=2*[137750:0.5:141084.5];
% InputR_DepTest=2*[141085:0.5:142084.5];
% InputR_DepPost=2*[142085:0.5:142749.5];
% Delta=2*[147750:0.5:147759.5];
% Ramp=2*[152760:0.5:155259.5];
% Discharge=[];
% for b=0:10
%     DischargeStart=160260+b*6000;
%     DischargeStop=160259+b*6000+1000;
%     Discharge(b+1,:)=2*[DischargeStart:0.5:DischargeStop+0.5];
% end

%File indices BY ERV fixed for AP
Baseline=[1:5000];
AP_Drop=[567020:572020];
AP_Waveform=[582020:582520];
% ERV IV
% IV=repmat((10000:20000:210000)', 1, 11000)+repmat(0:1:10999, 11,1);
% IV(6,:)=[]; %Eliminates the "no current" section
% AdamP IV
IV=[];
for a=0:10
    IVstart=5000+a*10000;
    IVstop=10000+a*10000;
    IV(a+1,:)=2*[IVstart:0.5:IVstop];
end
InputR_Hyp=[230005:240005];
InputR_HypTest=[236675:238675];
InputR_Dep=2*[250005:260005];
InputR_DepTest=2*[256674:258675];
Delta=[270030:272000];
PulseBeforeRamp=[295025:295525];
Ramp=[320000:325030];
% ERV Discharge
% Discharge=[[335025 347025:15000:482025]' [337025 352025:15000:487025]'];
% AdamP Discharge
Discharge=[];
for b=0:9
        DischargeStart=347025+b*15000;
        DischargeStop=352025+b*15000;
    Discharge(b+1,:)=DischargeStart:1:DischargeStop;
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find resting membrane potential if abs(holding current) < 5 pA
% if abs(mean(data(1:9999,2)))>5
%     Ephys.RestingPotential=[];
% else
%     Ephys.RestingPotential=mean(data(1:9999,1));
% end

%Baseline or average resting potential BY ERV
try
    Ephys.Baseline=mean(data(Baseline)); %THIS IS AN OUTPUT
    Ephys.BaselineShift=[mean(data(Baseline(1:5000))) mean(data(Baseline(5001:end)))]; %THIS IS AN OUTPUT
catch
    Ephys.Baseline=str2double('NaN');
    Ephys.BaselineShift=str2double('NaN');
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AP_Drop parameters
% Find APs and get their amplitudes
AP_Drop_ApIds=findaps2(data(AP_Drop,1));
AP_Drop_Amps=data(AP_Drop(AP_Drop_ApIds));
Ephys.InitAPDrop=AP_Drop_Amps(1)-AP_Drop_Amps(2);
Ephys.AP1ToSteadyDrop=AP_Drop_Amps(1)-AP_Drop_Amps(end);
Ephys.AP2ToSteadyDrop=AP_Drop_Amps(2)-AP_Drop_Amps(end);

% Find the change in amplitudes and get the absolute maximum change
AP_Drop_ChngAmps=diff(AP_Drop_Amps);
AP_Drop_MaxChngAmpIDX=find(max(abs(AP_Drop_ChngAmps)));
Ephys.MaxRateAPChange=AP_Drop_ChngAmps(AP_Drop_MaxChngAmpIDX);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AP_Waveform parameters
% Find when APs start (onset) and stop based on a change > 50 mV/ms
AP_Waveform_Onsets=find(diff(data(AP_Waveform,1))>5);
pts2kill=[];
for i=2:length(AP_Waveform_Onsets);    
    if AP_Waveform_Onsets(i)==AP_Waveform_Onsets(i-1)+1;
        pts2kill=[pts2kill i];
    end
end
AP_Waveform_Onsets(pts2kill)=[];
AP_Waveform_Onset_Voltages=data(AP_Waveform(AP_Waveform_Onsets));
AP_Waveform_Stops=[];
for i=1:length(AP_Waveform_Onsets);
    AP_Waveform_Stops=[AP_Waveform_Stops AP_Waveform_Onsets(i)+find(data(AP_Waveform(AP_Waveform_Onsets(i)):end,1)<AP_Waveform_Onset_Voltages(i),1)];
end

% Calculate Amplitude, Half Amplitude, and fAHP for AP 1 and 2
AP_Waveform1=data(AP_Waveform(AP_Waveform_Onsets(1):AP_Waveform_Stops(1)),1);
AP_Waveform2=data(AP_Waveform(AP_Waveform_Onsets(2):AP_Waveform_Stops(2)),1);
[AP_Waveform1_Amp,AP_Waveform1_AmpIdx]=max(AP_Waveform1);
[AP_Waveform2_Amp,AP_Waveform2_AmpIdx]=max(AP_Waveform2);
AP_Waveform1_HalfAmp=(AP_Waveform1_Amp-AP_Waveform_Onset_Voltages(1))/2;
AP_Waveform2_HalfAmp=(AP_Waveform2_Amp-AP_Waveform_Onset_Voltages(2))/2;
AP_Waveform1_HalfAmpOnIdx=find(AP_Waveform1(1:AP_Waveform1_AmpIdx)>AP_Waveform1(1)+AP_Waveform1_HalfAmp,1);
AP_Waveform2_HalfAmpOnIdx=find(AP_Waveform2(1:AP_Waveform2_AmpIdx)>AP_Waveform2(1)+AP_Waveform2_HalfAmp,1);
AP_Waveform1_HalfAmpOffIdx=AP_Waveform1_AmpIdx + find(AP_Waveform1(AP_Waveform1_AmpIdx:end)<AP_Waveform1(1)+AP_Waveform1_HalfAmp,1);
AP_Waveform2_HalfAmpOffIdx=AP_Waveform2_AmpIdx + find(AP_Waveform2(AP_Waveform2_AmpIdx:end)<AP_Waveform2(1)+AP_Waveform2_HalfAmp,1);
AP_Waveform1_fAHP=min(data(AP_Waveform(AP_Waveform_Stops(1):AP_Waveform_Onsets(2))));
try
    AP_Waveform2_fAHP=min(data(AP_Waveform(AP_Waveform_Stops(2):AP_Waveform_Onsets(3))));
catch
    AP_Waveform2_fAHP=min(data(AP_Waveform(AP_Waveform_Stops(2):end)));
end

% Write out params for AP 1
Ephys.AP1Amp=AP_Waveform1_Amp;
Ephys.AP1Duration=(AP_Waveform_Stops(1)-AP_Waveform_Onsets(1))/10;
Ephys.AP1HalfWidth=(AP_Waveform1_HalfAmpOffIdx-AP_Waveform1_HalfAmpOnIdx)/10;
Ephys.AP1RiseTime=(AP_Waveform1_AmpIdx-1)/10;
Ephys.AP1FallTime=(1+AP_Waveform_Stops(1)-(AP_Waveform1_AmpIdx+AP_Waveform_Onsets(1)))/10;
Ephys.AP1RiseRate=Ephys.AP1Amp/Ephys.AP1RiseTime;
Ephys.AP1FallRate=Ephys.AP1Amp/Ephys.AP1FallTime;
Ephys.AP1fAHP=AP_Waveform_Onset_Voltages(1)-AP_Waveform1_fAHP;

% Write out params for AP 2
Ephys.AP2Amp=AP_Waveform2_Amp;
Ephys.AP2Duration=(AP_Waveform_Stops(2)-AP_Waveform_Onsets(2))/10;
Ephys.AP2HalfWidth=(AP_Waveform2_HalfAmpOffIdx-AP_Waveform2_HalfAmpOnIdx)/10;
Ephys.AP2RiseTime=(AP_Waveform2_AmpIdx-1)/10;
Ephys.AP2FallTime=(1+AP_Waveform_Stops(2)-(AP_Waveform2_AmpIdx+AP_Waveform_Onsets(2)))/10;
Ephys.AP2RiseRate=Ephys.AP2Amp/Ephys.AP2RiseTime;
Ephys.AP2FallRate=Ephys.AP2Amp/Ephys.AP2FallTime;
Ephys.AP2fAHP=AP_Waveform_Onset_Voltages(2)-AP_Waveform2_fAHP;

% Write out params for changes between AP 1 and AP 2
Ephys.AP12AmpPercChng=(Ephys.AP1Amp-Ephys.AP2Amp)/Ephys.AP1Amp*100;
Ephys.AP12DurationPercChng=(Ephys.AP1Duration-Ephys.AP2Duration)/Ephys.AP1Duration*100;
Ephys.AP12HalfWidthPercChng=(Ephys.AP1HalfWidth-Ephys.AP2HalfWidth)/Ephys.AP1HalfWidth*100;
Ephys.AP12RiseRatePercChng=(Ephys.AP1RiseRate-Ephys.AP2RiseRate)/Ephys.AP1RiseRate*100;
Ephys.AP12FallRatePercChng=(Ephys.AP1FallRate-Ephys.AP2FallRate)/Ephys.AP1FallRate*100;
Ephys.AP12fAHPPercChng=(Ephys.AP1fAHP-Ephys.AP2fAHP)/Ephys.AP1fAHP*100;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IV
% Voltage peaks are min/max depending on whether hyperpol/depolar injection
% Voltage steady state is last 100ms of injection
% Current is mean of the current injection
IV_Vpeak=zeros(1,11);
IV_Vsteady=zeros(1,11);
IV_I=zeros(1,11);

% Negative current injections
for i=1:5
    IV_Vpeak(i)=min(data(IV(i,:),1));
    IV_Vsteady(i)=mean(data(IV(i,9000:10000),1));
    IV_I(i)=mean(data(IV(i,5000:9000),2));
end

% Positive current injections
for i=7:11
    IV_Vpeak(i)=max(data(IV(i,:),1));
    IV_Vsteady(i)=mean(data(IV(i,9000:10000),1));
    IV_I(i)=mean(data(IV(i,5000:9000),2));
end

% Delete middle non-existent IV injection
IV_Vpeak(6)=[];
IV_Vsteady(6)=[];
IV_I(6)=[];

IV_Rpeak=(IV_Vpeak./IV_I(1))*1000;
IV_Rsteady=(IV_Vsteady./IV_I(1))*1000;
IV_Sag=IV_Vpeak-IV_Vsteady;

% Calculate input resistance from first pulse only
% Multiply by 1000 to convert to megaohms
Ephys.InputRPeak=1000*(IV_Vpeak(1)-Ephys.Baseline)/IV_I(1);
Ephys.InputRSteady=1000*(IV_Vsteady(1)-Ephys.Baseline)/IV_I(1);

Ephys.RectificationPeak=(IV_Rpeak(1)-abs(IV_Rpeak(10)))/IV_Rpeak(1);
Ephys.RectificationSteady=(IV_Rsteady(1)-abs(IV_Rsteady(10)))/IV_Rsteady(1);
Ephys.Sag=abs(max(abs(IV_Sag)));

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Delta
% Fit voltage response to a short current step to a single exponential 
% in order to get membrane time constant tau
% This calculation is probably faulty because the 'capacitor' (in this
% case, the membrane) will not have fully charged within such a short time!

% AdamP code
% DeltaV=data(Delta(:,1));
% DeltaT=(0.1:0.1:length(DeltaV)/10)';
% betaa=nlinfit(DeltaT,DeltaV-max(DeltaV),'acq_single_exp',[min(DeltaV) 1]);
% Ephys.DeltaTau=betaa(2);

% ERV code
DeltaV=data(Delta(1,:));
FlippedDeltaV=-DeltaV;
ToFit=FlippedDeltaV-min(FlippedDeltaV);
[yhat, Ephys.tau]=single_exp(ToFit,10000);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ramp
Ramp_APIdx=findaps2(data(Ramp,1));
Ramp_Thresh=data(Ramp(Ramp_APIdx(1)),2);
Ramp_AHP1=min(data(Ramp(Ramp_APIdx(1):Ramp_APIdx(2))));
Ramp_AP_DiffGreater5=find(diff(data(Ramp,1))>5);
Ramp_Onset1Voltage=data(Ramp(Ramp_AP_DiffGreater5));

Ephys.RampThresh=Ramp_Thresh;
Ephys.RampfAHP=Ramp_Onset1Voltage(1)-Ramp_AHP1;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sAHP/Discharge
% Get info about APs and AHP
Pattern_AllISIs=[];
for i=1:10
    Pattern_APIdx{i}=findaps2(data(Discharge(i,:),1));
    Pattern_AP1Delay(i)=(Pattern_APIdx{i}(1)-5000)/10;
    Pattern_AP2Delay(i)=(Pattern_APIdx{i}(2)-5000)/10;
    Pattern_ISIs{i}=diff(Pattern_APIdx{i});
    Pattern_MeanOfISIs123(i)=(mean(Pattern_ISIs{i}(1:3)))/10;
    Pattern_InitISIChange(i)=(Pattern_ISIs{i}(2)-Pattern_ISIs{i}(1))/10;
    Pattern_LastISIChange(i)=(Pattern_ISIs{i}(end)-Pattern_ISIs{i}(end-1))/10;
    DischargeStop=Discharge(i,end);
    [Pattern_AHP(i) Pattern_AHPIdx(i)]=min(data(DischargeStop:DischargeStop+10000,1));
    % If the last AP occurs after current injection, set the Offset Time
    % for the burst equal to that last AP's index
    if Pattern_APIdx{i}(end) > 7000
        Pattern_OffsetIdx(i)=Pattern_APIdx{i}(end);
    else
        Pattern_OffsetIdx(i)=5000;
    end
%     Pattern_100msAHP(i)=data(Discharge(i,Pattern_OffsetIdx(i)+1000),1);
%     Pattern_Baseline(i)=mean(data(Discharge(i,1:5000),1));
%     Pattern_InjCurrent(i)=mean(data(Discharge(i,5000:7000),2));
    Pattern_NumAPs(i)=length(Pattern_APIdx{i});
    Pattern_AllISIs=[Pattern_AllISIs Pattern_ISIs{i}/10];
    Pattern_MaxDerivMinusAveDerivISIs(i)=max(diff(Pattern_ISIs{i}))-mean(diff(Pattern_ISIs{i}));
end
% Pattern_sAHP=Pattern_Baseline-Pattern_AHP;
% Ephys.PostBurstMaxAHP=max(Pattern_sAHP);
% Ephys.PostBurst100msAHP=mean(Pattern_100msAHP-Pattern_Baseline);
% Ephys.PostBurstTimeToMaxAHP=mean(((Pattern_AHPIdx+7000)-Pattern_OffsetIdx)/10);
% Ephys.NumSpikesPerPicoAmp=mean(Pattern_NumAPs)/mean(Pattern_InjCurrent);
Ephys.AveDelayToFirstSpike=mean(Pattern_AP1Delay);
Ephys.StdDelayToFirstSpike=std(Pattern_AP1Delay);
Ephys.AveDelayToSecondSpike=mean(Pattern_AP2Delay);
Ephys.StdDelayToSecondSpike=std(Pattern_AP2Delay);
Ephys.AveFirstThreeApISIs=mean(Pattern_MeanOfISIs123);
Ephys.SDFirstThreeApISIs=std(Pattern_MeanOfISIs123);
Ephys.AveInitFiringRateAccom=mean(Pattern_InitISIChange);
Ephys.AveSteadyFiringRateAccom=mean(Pattern_LastISIChange)-Ephys.AveInitFiringRateAccom;
Ephys.ISI_CV=std(Pattern_AllISIs)/mean(Pattern_AllISIs);
Ephys.ISIMedian=median(Pattern_AllISIs);
Ephys.aveMaxDerivMinusAveDerivISIs=(mean(Pattern_MaxDerivMinusAveDerivISIs))/10;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sub-Subfunction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function yhat=single_exp_fun(beta, t)
% f0=beta(1);
% tau=beta(2);
% yhat=f0.*exp(-t./tau);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Subfunction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [yhat, tau]=single_exp(data, samplerate)
single_exp_fun = @(beta,t) beta(1).*exp(-t./beta(2));
t=(0:1/samplerate:(length(data)/samplerate)-(1/samplerate)).*1000; %Time scale in milliseconds
FittedBetas=nlinfit(t, data,single_exp_fun, [min(data) 1]);
yhat=single_exp_fun(FittedBetas, t);
tau=FittedBetas(2);
