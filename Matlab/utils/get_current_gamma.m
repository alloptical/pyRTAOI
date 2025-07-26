function [this_gamma,this_xi,alpha,L] = get_current_gamma(this_data,hmm,alpha_old)
% adapted from hsinference in hmm-mar package for online analysis
% only works for order = 0; covtype = 'uniquefull' 
% ZZ 2019

K = hmm.K; % number of states
P = hmm.P; % the (K X K) state transition probability matrices
Pi = hmm.Pi; % the (K X 1) initial state probabilities.

%% ahared params - move to input
%these params should be saved in cache as they don't change with state or data
S = hmm.train.S==1;
regressed = sum(S,1)>0;
ltpi = sum(regressed)/2 * log(2*pi);
ldetWishB=0.5*logdet(hmm.Omega.Gam_rate(regressed,regressed));
PsiWish_alphasum=0;
for n=1:sum(regressed)
    PsiWish_alphasum=PsiWish_alphasum+psi(hmm.Omega.Gam_shape/2+0.5-n/2);
end
PsiWish_alphasum=PsiWish_alphasum*0.5;
C = hmm.Omega.Gam_shape * hmm.Omega.Gam_irate;
NormWishtraces = zeros(1,K);
meands = {};
for k = 1:K
    NormWishtraces(k) = 0.5 * sum(sum(C .* hmm.state(k).W.S_W));
    meands{k}= hmm.state(k).W.Mu_W(:,regressed);
end

%% get observation likelihood for each state
% adapted from obslike for current frame with default hmm params
L = zeros(1,hmm.K);
for k = 1:hmm.K    % loop through states
    d = this_data(regressed) - meands{k};
    dist = - 0.5*d*C*d';
    L(k)= - ltpi - ldetWishB + PsiWish_alphasum + dist - NormWishtraces(k);
end
L = exp(L);
L(L<realmin) = realmin;
%% get inference using forward algorithm
if isempty(alpha_old)
    alpha = Pi.*L;
else
    alpha = (alpha_old.*P).*L;
end
scale = sum(alpha);	
alpha = alpha./scale;
t = P.*( alpha');
this_gamma = alpha;
this_xi =  t(:)'/sum(t(:));
end

