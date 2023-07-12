function [z] = resVarSASS(y,lambdas,type)
% denoising using SASS with method of residual variances for lambda tuning 
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
fc = 0.05;     % fc : cut-off frequency (cycles/sample) (0 < fc < 0.5);
K = 2;          % K : order of sparse derivative

Errs = zeros([1, length(lambdas)]); % store Error values
% Compute residual variance for each parameter value
for i = 1:length(lambdas)
    lam = lambdas(i);
    [f, ~, ~, ~, ~, ~] = sass_L1(y, d, fc, K, lam);
    Errs(:,i) = (sum((y-f).^2))/length(y);
end

%Choose lambda
[~,i] = min(abs(Errs - V));
[z, ~, ~, ~, ~, ~] = sass_L1(y, d, fc, K, lambdas(i));
end