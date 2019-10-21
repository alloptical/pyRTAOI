function Ephys=packerwave(data, wave)


%Parameters
acqrate=10; %Acquisition rate in KHz
commandgain=400; %External command sensitivity from the multiclamp amplifier 400pA/V (default) or 2000pA/V
AP_threshold=-20; %AP threshold is -20mV
AP_space=5; %space between APs is 5ms
AP_speed=70; %Speed threshold to find onset in mV/ms
IV_steady=100; %Miliseconds to average to calculate steady state potential during the IV
IV_extra=100; %Miliseconds after the current injections for the IV finish

%File indeces
Baseline=[1:5000];
AP_Drop=[567020 572020];
AP_Waveform=[582020 582520];
IV=repmat((10000:20000:210000)', 1, 11000)+repmat(0:1:10999, 11,1);
IV(6,:)=[]; %Eliminates the "no current" section
InputR_Hyp=[230005 240005];1
InputR_HypTest=[236675 238675];
InputR_Dep=2*[250005 260005];
InputR_DepTest=2*[256674 258675];
Delta=[270002 270026];
PulseBeforeRamp=[295025 295525];
Ramp=[320000 325030];
Discharge=[[335025 347025:15000:482025]' [337025 352025:15000:487025]'];

%Baseline or average resting potential
try
    Ephys.Baseline=mean(data(Baseline)); %THIS IS AN OUTPUT
    Ephys.BaselineShift=[mean(data(Baseline(1:5000))) mean(data(Baseline(5001:end)))]; %THIS IS AN OUTPUT
catch
    Ephys.Baseline=str2double('NaN');
    Ephys.BaselineShift=str2double('NaN');
end

%AP detection for "drop" calculations
datatemp=data(AP_Drop(1):AP_Drop(2));
a=find(datatemp>AP_threshold); %AP threshold
[index, b]=find(diff(a)>(AP_space*acqrate)); %space between APs
index=[0; index; length(a)];
Ephys.APNumber=length(index)-1; %THIS IS AN OUTPUT
Ephys.APAmplitude=zeros(length(index)-1,1);
for i=1:length(index)-1
    Ephys.APAmplitude(i)=max(datatemp(a(index(i)+1):a(index(i+1)))); %THIS IS AN OUTPUT
end
Ephys.APAmplitude=Ephys.APAmplitude-Ephys.Baseline;
Ephys.InitAPDrop=Ephys.APAmplitude(1)-Ephys.APAmplitude(2);
Ephys.AP1ToSteadyDrop=Ephys.APAmplitude(1)-Ephys.APAmplitude(end);
Ephys.AP2ToSteadyDrop=Ephys.APAmplitude(2)-Ephys.APAmplitude(end);
Ephys.MaxRateAPChange=max(abs(diff(Ephys.APAmplitude)));


%AP parameters for the first AP triggered by the last pulse (waveform in Adam's nomenclature)
datatemp=data(AP_Waveform(1):AP_Waveform(2));
a=find(diff(datatemp)>AP_speed/acqrate); %Start of AP
b=find(datatemp(a(1):end)<datatemp(a(1))); %End of AP
[i, index]=max(datatemp(a(1)+b(1)-1:end));
i=min(datatemp(a(1)+b(1)-1:a(1)+b(1)-1+index));
datatemp=(datatemp(a(1):a(1)+b(1)-1));
[Ephys.AP1Amp, c]=max(datatemp); %THIS IS AN OUTPUT
Ephys.AP1AmpAbs=Ephys.AP1Amp-Ephys.Baseline; %THIS IS AN OUTPUT
Ephys.AP1Duration=(b(1)-1)/acqrate; %THIS IS AN OUTPUT
d=(Ephys.AP1Amp-datatemp(1))/2;
d=Ephys.AP1Amp-d;
Ephys.AP1HalfWidth=(interp1(datatemp(c:end),c:length(datatemp),d)-interp1(datatemp(1:c),1:c,d))/acqrate; %THIS IS AN OUTPUT
Ephys.AP1RiseTime=(c-1)/acqrate; %THIS IS AN OUTPUT
Ephys.AP1FallTime=(length(datatemp)-c)/acqrate; %THIS IS AN OUTPUT
Ephys.AP1RiseRate=Ephys.AP1Amp/Ephys.AP1RiseTime; %THIS IS AN OUTPUT
Ephys.AP1FallRate=Ephys.AP1Amp/Ephys.AP1FallTime; %THIS IS AN OUTPUT
Ephys.AP1fAHP=datatemp(1)-i; %THIS IS AN OUTPUT


%Average AP parameters from first AP triggered by the pulse before the ramp and the 11 pulses after the ramp
Ephys.AP1AmpAv=[];
Ephys.AP1AmpAbsAv=[];
Ephys.AP1DurationAv=[];
Ephys.AP1HalfWidthAv=[];
Ephys.AP1RiseTimeAv=[];
Ephys.AP1FallTimeAv=[];
Ephys.AP1RiseRateAv=[];
Ephys.AP1FallRateAv=[];
Ephys.AP1fAHPAv=[];
ranges=[PulseBeforeRamp; Discharge];
for i=1:12
    datatemp=data(ranges(i,1):ranges(i,2));
    a=find(diff(datatemp)>AP_speed/acqrate); %Start of AP
    b=find(datatemp(a(1):end)<datatemp(a(1))); %End of AP
    [i, index]=max(datatemp(a(1)+b(1)-1:end));
    i=min(datatemp(a(1)+b(1)-1:a(1)+b(1)-1+index));
    datatemp=(datatemp(a(1):a(1)+b(1)-1));
    [e, c]=max(datatemp);
    Ephys.AP1AmpAv=[Ephys.AP1AmpAv e];
    Ephys.AP1AmpAbsAv=[Ephys.AP1AmpAbsAv Ephys.AP1AmpAv(end)-Ephys.Baseline];
    Ephys.AP1DurationAv=[Ephys.AP1DurationAv (b(1)-1)/acqrate];
    d=(Ephys.AP1AmpAv(end)-datatemp(1))/2;
    d=Ephys.AP1AmpAv(end)-d;
    Ephys.AP1HalfWidthAv=[Ephys.AP1HalfWidthAv (interp1(datatemp(c:end),c:length(datatemp),d)-interp1(datatemp(1:c),1:c,d))/acqrate];
    Ephys.AP1RiseTimeAv=[Ephys.AP1RiseTimeAv (c-1)/acqrate];
    Ephys.AP1FallTimeAv=[Ephys.AP1FallTimeAv (length(datatemp)-c)/acqrate];
    Ephys.AP1RiseRateAv=[Ephys.AP1RiseRateAv Ephys.AP1AmpAv(end)/Ephys.AP1RiseTimeAv(end)];
    Ephys.AP1FallRateAv=[Ephys.AP1FallRateAv Ephys.AP1AmpAv(end)/Ephys.AP1FallTimeAv(end)];
    Ephys.AP1fAHPAv=[Ephys.AP1fAHPAv datatemp(1)-i];
end
Ephys.AP1AmpAv=mean(Ephys.AP1AmpAv); %THIS IS AN OUTPUT
Ephys.AP1AmpAbsAv=mean(Ephys.AP1AmpAbsAv); %THIS IS AN OUTPUT
Ephys.AP1DurationAv=mean(Ephys.AP1DurationAv); %THIS IS AN OUTPUT
Ephys.AP1HalfWidthAv=mean(Ephys.AP1HalfWidthAv); %THIS IS AN OUTPUT
Ephys.AP1RiseTimeAv=mean(Ephys.AP1RiseTimeAv); %THIS IS AN OUTPUT
Ephys.AP1FallTimeAv=mean(Ephys.AP1FallTimeAv); %THIS IS AN OUTPUT
Ephys.AP1RiseRateAv=mean(Ephys.AP1RiseRateAv); %THIS IS AN OUTPUT
Ephys.AP1FallRateAv=mean(Ephys.AP1FallRateAv); %THIS IS AN OUTPUT
Ephys.AP1fAHPAv=mean(Ephys.AP1fAHPAv); %THIS IS AN OUTPUT


%IV
datatemp=data(IV)';
IV_Vpeak=[min(datatemp(:,1:5)) max(datatemp(:,6:end))];
a=find(IV_Vpeak>-20);
IV_Vpeak=IV_Vpeak-Ephys.Baseline;
IV_Vsteady=mean(datatemp(end-(IV_extra+IV_steady)*acqrate:end-IV_extra*acqrate, :))-Ephys.Baseline;
datatemp=wave(IV)';
IV_I=mean(datatemp(10*acqrate:end-(IV_extra+10)*acqrate, :))*commandgain; %Gets rid of the first and las 10ms just in case
disp(['Pulse(s) ', num2str(a), ' contain APs']);
if ~isempty(a)
    IV_Vpeak(a)=[];
    IV_Vsteady(a)=[];
    IV_I(a)=[];
end
Ephys.InputRPeakHyp=polyfit(IV_I(1:5),IV_Vpeak(1:5),1); %THIS IS AN OUTPUT
Ephys.InputRPeakHyp=Ephys.InputRPeakHyp(1)*1000; %Multiplied by 1000 to convert GOhms in MOhms %THIS IS AN OUTPUT
Ephys.InputRPeakDep=polyfit(IV_I(6:end),IV_Vpeak(6:end),1); %THIS IS AN OUTPUT
Ephys.InputRPeakDep=Ephys.InputRPeakDep(1)*1000; %Multiplied by 1000 to convert GOhms in MOhms %THIS IS AN OUTPUT
Ephys.RectificationPeak=Ephys.InputRPeakHyp/Ephys.InputRPeakDep; %THIS IS AN OUTPUT
Ephys.InputRSteadyHyp=polyfit(IV_I(1:5),IV_Vsteady(1:5),1); %THIS IS AN OUTPUT
Ephys.InputRSteadyHyp=Ephys.InputRSteadyHyp(1)*1000; %Multiplied by 1000 to convert GOhms in MOhms %THIS IS AN OUTPUT
Ephys.InputRSteadyDep=polyfit(IV_I(6:end),IV_Vsteady(6:end),1); %THIS IS AN OUTPUT
Ephys.InputRSteadyDep=Ephys.InputRSteadyDep(1)*1000; %Multiplied by 1000 to convert GOhms in MOhms %THIS IS AN OUTPUT
Ephys.RectificationSteady=Ephys.InputRSteadyHyp/Ephys.InputRSteadyDep; %THIS IS AN OUTPUT
Ephys.Sag=abs(IV_Vpeak(1))-abs(IV_Vsteady(1)); %THIS IS AN OUTPUT


%Ramp
datatemp=data(Ramp(1):Ramp(2));
a=find(diff(datatemp)>AP_speed/acqrate); %Start of AP
b=find(diff(a)>1);
datatemp=datatemp(a(1):a(b(1)+1));
Ephys.RampThresh=datatemp(1);
Ephys.RampfAHP=datatemp(1)-min(datatemp);

%Membrane time constant using the first part of the pulse (the falling part)
datatemp=data(Delta(1):Delta(2));
[i, index]=min(datatemp);
datatemp=datatemp(1:index);
% [yhat, Ephys.tau]=single_exp(datatemp, acqrate*1000);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Subfunction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [yhat, tau]=single_exp(data, samplerate)
data=data-min(data);
single_exp_fun = @(beta,t) beta(1).*exp(-t./beta(2));
t=(0:1/samplerate:(length(data)/samplerate)-(1/samplerate)).*1000; %Time scale in milliseconds
FittedBetas=nlinfit(t, data,single_exp_fun, [min(data) 1]);
yhat=single_exp_fun(FittedBetas, t);
tau=FittedBetas(2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sub-Subfunction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function yhat=single_exp_fun(beta, t)
% f0=beta(1);
% tau=beta(2);
% yhat=f0.*exp(-t./tau);


