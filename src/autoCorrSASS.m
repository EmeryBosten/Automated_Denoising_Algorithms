function [z] = autoCorrSASS(y,lambdas,type)
% denoising using SASS with method of autocorrelation for lambda tuning 
% Input
% y : noisy signal
% lambdas : set of regularization parameters
% Output 
% z : denoised signal

% compute autocorrelation value 
e = diff(y);
if strcmp(type,'mean')
    p = autocorrelation(e,'mean');
elseif strcmp(type,'median')
    p = autocorrelation(e,'percentile');
end

d = 2;          % d : filter order parameter (d = 1, 2, or 3)
fc = 0.05;      % fc : cut-off frequency (cycles/sample) (0 < fc < 0.5);
K = 2;          % K : order of sparse derivative

Ps = zeros([1, length(lambdas)]);
for i = 1:length(lambdas)
    lam = lambdas(i);
    [z, ~, ~, ~, ~, ~] = sass_L1(y, d, fc, K, lam);
    Ps(:,i) = autocorrelation(y-z,'mean');
end

%Choose lambda
[minval,i] = min(abs((Ps - p)));
[z, ~, ~, ~, ~, ~] = sass_L1(y, d, fc, K, lambdas(i));
end
