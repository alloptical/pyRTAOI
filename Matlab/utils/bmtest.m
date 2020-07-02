function [b] = bmtest(x)
% copied from Freeman and Dale
% Assessing bimodality to detect the presence of a dual cognitive process
% bimodality test
% a bimodal
% distribution will have very low kurtosis, an asymmetric
% character, or both; all of these conditions increase BC. The
% values range from 0 and 1, with those exceeding .555 (the
% value representing a uniform distribution) suggesting bimodality (SAS Institute, 1989)

%m3 = skew
%m4 = kurt
%n = data size
m3 = skewness(x);
m4 = kurtosis(x);
n = length(x);
b=(m3^2+1) / (m4 + 3 * ( (n-1)^2 / ((n-2)*(n-3)) )); 
end

