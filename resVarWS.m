function [z] = resVarWS(y,lambdas,type)
% denoising using WS with method of residual variances for lambda tuning 
% Input
% y : noisy signal
% lambdas : set of regularization parameters
% Output 
% z : denoised signal

%Compute estimate of redidual variance of noisy signal y by modified one 
% nearest neighbor estimator
x=linspace(1,length(y),length(y))';
[V] = modoneNNestimator(y,x,type);

d = 2;          % d : filter order parameter (d = 1, 2, or 3)
Errs = zeros([1, length(lambdas)]);
for i = 1:length(lambdas)
    lam = lambdas(i);
    [z] = WSsparse(y', lam, d);
    Errs(:,i) = (sum((y-z).^2))/length(y);
end

%Choose lambda
[~,i] = min(abs(Errs - V));
[z] = WSsparse(y', lambdas(i),d);
end