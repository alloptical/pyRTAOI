function yhat=acq_single_exp(beta, t)

% Function used to calculate the decay time constant
% of currents using the nlnfit function
% yhat are the fitted values
% beta is a vector with two initial guesses of the values: f0 and tau
% t is the time vector in msec.

f0=beta(1);
tau=beta(2);
yhat=f0.*exp(-t./tau);