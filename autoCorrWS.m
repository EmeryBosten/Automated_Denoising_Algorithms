function [z] = autoCorrWS(y,lambdas,type)
% denoising using WS with method of autocorrelation for lambda tuning 
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
Ps = zeros([1, length(lambdas)]);
for i = 1:length(lambdas)
    lam = lambdas(i);
    [z] = WSsparse(y', lam, d);
    Ps(:,i) = autocorrelation(y-z,type);
end

%Choose lambda
[minval,i] = min(abs((Ps - p)));
[z] = WSsparse(y', lambdas(i),d);
end