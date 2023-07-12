function [z,b] = autoCorrBEADS(y,alphas,type)
% denoising and baseline correction using BEADS with method of autocorrelation for alpha tuning 
% Input
% y : noisy signal
% lambdas : set of relation parameters
% Output 
% z : denoised signal
% b : baseline

% compute autocorrelation value 
e = diff(y);
if strcmp(type,'mean')
    p = autocorrelation(e,'mean');
elseif strcmp(type,'median')
    p = autocorrelation(e,'percentile');
end

fc = 0.006;             % fc : cut-off frequency (cycles/sample)
d = 1;                  % d : filter order parameter (d = small positive integer)
Nit = 50;
EPS0 = 1e-5;            % for x itself
EPS1 = 1e-5;            % for derivatives
r = 1;                  % to set the asymmetric ratio
pen = 'L1_v1';          % to choose penalty function

Ps = zeros([1, length(alphas)]);
for i = 1:length(alphas)
    alpha = alphas(i);
    [z,b,~,~,~] = beads(y, d, fc, r, alpha, pen, Nit, EPS0, EPS1);
    Ps(:,i) = autocorrelation(y-(z+b),type);
end

%Choose lambda
[minval,i] = min(abs((Ps - p)));
[z,b,~,~,~] = beads(y, d, fc, r, alphas(i), pen, Nit, EPS0, EPS1);
end