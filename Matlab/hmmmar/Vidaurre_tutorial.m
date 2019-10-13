% This tutorial requires the installation of two toolboxes
%
% OSL: for preprocessing of MEG data (in this case used for visualisation
% purposes)
% HMM-MAR: containing the code of the HMM (and, regardless of what the name
% suggests, implementing other models apart from the MAR)
%
% The script install.sh installs both
%
% Note that a tar.gz with OSL is provided, but can also be downloaded from
% https://ohba-analysis.github.io/osl-docs/pages/overview/download.html
% 
% If the code complains, then FSL must be installed as well. 
% For this, follow this link: https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FslInstallation
%
% The total running time of this tutorial (in my laptop) is around 10min,
% but I have included a file with the precomputed outputs in case it runs
% too slow. 

%% Prepare data

cd /Users/dvidaurre/Work/Matlab/Course_NY % Set it to your path
cd C:\Users\chris\Desktop\neural_data_analysis_NYC2019\Wanda code and data\Diego Vidaurre
softwaredir = pwd; %'./'; % adjust to your path

outdir = 'out/';
subj_str = 'subj_2.mat';

% Check the atlas file: which regions to use? 
atlasfile = './parcellation/fmri_d100_parcellation_with_PCC_tighterMay15_v2_8mm.nii.gz';
addpath(genpath([softwaredir 'osl.tar/osl/osl-core'])) % OSL toolbox
addpath(genpath(pwd)) % OSL toolbox
osl_startup
p = parcellation(atlasfile);
osleyes(p);

addpath(genpath([softwaredir 'HMM-MAR-master'])) % HMM toolbox
addpath('scripts')

channels = [1 2 3 4 13 14 16 17 26 27 34 35]; % which channels to run the HMM on
conditions = {'Famous face','Unfamiliar face','Scrambled face'};

%% Run unsupervised HMM

datadir = 'data/'; % Here the data (X) is organised as continuous data (time points by channels)

options = struct(); % prepare HMM options
options.standardise = 1;     % standardize the data *only* in supervised analysis
options.Fs = 250;            % sampling frequency
options.K = 6;  	         % The number of states to infer
options.order = 0; 	         % The lag used, this is only relevant when using MAR observations
options.zeromean = 1; 	     % We do not want to model the mean, so zeromean is set on
options.covtype = 'full';    % We want to model the full covariance matrix
options.embeddedlags = -7:7; % 15 lags are used from -7 to 7
options.pca = length(channels)*2;   % The PCA dimensionality reduction is 2 times the number of ROIs 
%%% because the covar matrix bins x channels is high
options.initrep = 1; % to make it quicker - leave by default otherwise
options.initcyc = 1; % to make it quicker - leave by default otherwise 
options.cyc = 20; % to make it quicker - leave by default otherwise
options.verbose = 1; 
options_copy = rmfield(options,'embeddedlags');
window = 2;

load([datadir subj_str]); % 50 min of MEG data for one subject, 38 ROIs
tic
[hmm,Gamma_hmm] = hmmmar(X(:,channels),T,options);
toc % 4 min
evokedGamma = cell(3,1);
evokedField = cell(3,1);
for c = 1:3 % three types of stimulus: 'Famous','Unfamiliar','Scrambled' faces
    stim = stimulus == c;
    evokedGamma{c} = evokedStateProbability(stim,T,Gamma_hmm,window,options);
    evokedField{c} = evokedStateProbability(stim,T,zscore(X(:,channels)),window,options_copy);
end
save([outdir subj_str],'hmm','Gamma_hmm','evokedGamma','evokedField')

%% Obtain spectral information of the states

spec_options = struct();
spec_options.fpass = [1 40]; % which band we care about
spec_options.win = 256; % window of the multitaper estimation 
spec_options.embeddedlags = -7:7; % set to same than in the HMM run
spec_options.Fs = 250; % sampling frequency

tic
spectra = hmmspectramt(X(:,channels),T,Gamma_hmm,spec_options);
toc % 1 minute
save([outdir subj_str],'spectra','-append')

%% Show evoked field response

t = -250:250;

figure(2)
load([outdir subj_str],'evokedField')
ch = 10; % pick one channel  

for c = 1:3
    subplot(1,3,c)
    plot(t/250,evokedField{c}(:,ch),'LineWidth',2); xlim([-1 1])
    xlabel('Time')
    title(conditions{c})
end


%% Show evoked state response

t = -250:250;

figure(3)
load([outdir subj_str],'evokedGamma')

for c = 1:3
    subplot(1,3,c)
    plot(t/250,evokedGamma{c},'LineWidth',2); xlim([-1 1])
    xlabel('Time')
    title(conditions{c})
    legend('state 1','state 2','state 3','state 4','state 5','state 6','Location','South')
end



%% Show spectra

figure(4);clf(4)
subset_channels = 1:4; 

cm = [0 0.4470 0.7410;
    0.8500 0.3250 0.0980;
    0.9290 0.6940 0.1250;
    0.4940 0.1840 0.5560;
    0.4660 0.6740 0.1880;
    0.3010 0.7450 0.9330];

load([outdir subj_str],'spectra')
plot_hmmspectra (spectra,[],[],4, cm, subset_channels);

%% Testing on state occupancy differences between conditions

window = [-0.2 0.5] ; % in seconds
samp_rate = 250;

[ trial_data, design_mat, time ] = load_epoch_results( Gamma_hmm, stimulus, window, samp_rate);
% trial_data has nsamples by states by ntrials
% design_mat has ntrials by conditions (famous, unfamiliar, scrambled)
nsamples = size(trial_data,1); 
K = size(Gamma_hmm,2); % no. of HMM states

ncontrasts = 2;             
contrasts = [   [0 1 -1 0];  ... % famous vs unfamiliar
                [0 1 1 -2] ...  % face vs non-face 
                ];

tic
pvals = zeros(ncontrasts,nsamples);
for t = 1:nsamples
    data = squeeze(trial_data(t,:,:))';
    good_inds = ~isnan(sum(data,2));
    for jj = 1:2
        Yin = design_mat(good_inds,:) * contrasts(jj,:)';
        Xin = data(good_inds,:);
        pvals(jj,t) = permregress(Xin,Yin,5000);
    end
    if mod(t,10)==0
        disp(['Time point ' num2str(t)])
    end
end
toc % 1min30s
save([outdir subj_str],'pvals','-append')


%% Show results for testing

load([outdir subj_str],'pvals')
figure(5);clf(5)
plot(time,pvals,'LineWidth',3)
hold on; plot(time,0.05*ones(1,length(time)),'k'); hold off
legend('Contrast 1','Contrast 2','significance')
xlabel('Time'); ylabel('P-value')

%% Show states in osleyes

net_mean = zeros(39,hmm.train.K);
for k = 1:length(spectra.state)
    net_mean(channels,k) = diag(squeeze(mean(abs(spectra.state(k).psd),1)));
end
p.osleyes(net_mean);

%% Run standard decoding on the data

datadir = 'data_epoched/'; % Here the data (X) is organised as epoched data (time points by trials by channels)
load([datadir subj_str])

% famous vs unfamiliar
get_out = Y(:,3)==1; % those that are not faces
Y2 = Y(~get_out,:) * contrasts(1,2:end)' ;
X2 = X(:,~get_out,:);
T2 = T(~get_out);
acc_famous_vs_familiar = standard_decoding(X2,Y2,T2);
baseline_famous_vs_familiar = min(mean(Y2==1),mean(Y2==-1));

% faces vs nonfaces
Y2 = Y * contrasts(2,2:end)' ;
acc_faces_vs_nonfaces = standard_decoding(X,Y2,T);
baseline_faces_vs_nonfaces = min(mean(Y2==1),mean(Y2==-2));

save([outdir subj_str],'acc_famous_vs_familiar','acc_faces_vs_nonfaces',...
    'baseline_famous_vs_familiar','baseline_faces_vs_nonfaces','-append')

%% Show cross-validated accuracy

figure(7);clf(7)

t = (1:751) - 250;

load([outdir subj_str],'acc_famous_vs_familiar','acc_faces_vs_nonfaces')
cv_acc = [acc_famous_vs_familiar acc_faces_vs_nonfaces];
plot(t/250,cv_acc,'LineWidth',2); 
hold on
plot(t,ones(size(t))*baseline_famous_vs_familiar,'b')
plot(t,ones(size(t))*baseline_faces_vs_nonfaces,'r')
hold off
xlim([t(1) t(end)]/250)
xlabel('Time')


%% Run supervised HMM

datadir = 'data_epoched/'; % Here the data (X) is organised as epoched data (time points by trials by channels)
load([datadir subj_str])

options = struct(); % prepare TUDA options
options.K = 4;
options.DirichletDiag = 100;
options.Fs = 250;
options.initrep = 1; % to make it quicker - leave by default otherwise
options.initcyc = 1; % to make it quicker - leave by default otherwise 
options.cyc = 20; % to make it quicker - leave by default otherwise

t = 250:650; % time points to use, based on unsupervised analysis and standard decoding

% first contrast
get_out = Y(:,3)==1; % those that are not faces
Y2 = Y(~get_out,:) * contrasts(1,2:end)' ;
X2 = X(t,~get_out,channels); % X is here (time points by trials by channels)
X2 = reshape(X2,size(X2,1)*size(X2,2), size(X2,3)); % concatenate time points and trials
T2 = T(~get_out);
T2 = length(t) * ones(size(T2)); % length of each trial in time points
tic
[tuda_famous_vs_familiar,Gamma_tuda_famous_vs_familiar] = tudatrain(X2,Y2,T2,options);
toc % 45s

% second contrast
Y2 = Y * contrasts(2,2:end)' ;
X2 = X(t,:,channels); % X is here (time points by trials by channels)
X2 = reshape(X2,size(X2,1)*size(X2,2), size(X2,3)); % concatenate time points and trials
T2 = length(t) * ones(size(T)); % length of each trial in time points
tic
[tuda_faces_vs_nonfaces,Gamma_tuda_faces_vs_nonfaces] = tudatrain(X2,Y2,T2,options);
toc % 1 min

save([outdir subj_str],'tuda_faces_vs_nonfaces','Gamma_tuda_faces_vs_nonfaces',...
    'tuda_famous_vs_familiar','Gamma_tuda_famous_vs_familiar','-append')


%% Show average progression of decoders

t = (250:650) - 250;

figure(8);clf(8)

% First contrast
subplot(211)
load([outdir subj_str],'Gamma_tuda_famous_vs_familiar','tuda_famous_vs_familiar')
Gamma_tuda = reshape(Gamma_tuda_famous_vs_familiar,[length(t) ...
    size(Gamma_tuda_famous_vs_familiar,1)/length(t) tuda_famous_vs_familiar.train.K]);
mGamma = squeeze(mean(Gamma_tuda,2));
plot(t/250,mGamma,'LineWidth',2); xlim([t(1) t(end)]/250)
xlabel('Time')
legend('state 1','state 2','state 3','state 4')

% Second contrast
subplot(212)
load([outdir subj_str],'tuda_faces_vs_nonfaces','Gamma_tuda_faces_vs_nonfaces')
Gamma_tuda = reshape(Gamma_tuda_faces_vs_nonfaces,[length(t) ...
    size(Gamma_tuda_faces_vs_nonfaces,1)/length(t) tuda_faces_vs_nonfaces.train.K]);
mGamma = squeeze(mean(Gamma_tuda,2));
plot(t/250,mGamma,'LineWidth',2); xlim([t(1) t(end)]/250)
xlabel('Time')
legend('state 1','state 2','state 3','state 4')


%% Show decoders in osleyes

beta = tudabeta(tuda_faces_vs_nonfaces);
beta = squeeze(sum(beta.^2,2)); % sum across conditions
beta_wholebrain = zeros(39,tuda_faces_vs_nonfaces.train.K); % no. parcels y states
beta_wholebrain(channels,:) = beta;
p.osleyes(beta_wholebrain);
