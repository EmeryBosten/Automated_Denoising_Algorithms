function [z,b] = resVarBEADS(y,alphas,type)
% denoising and baseline correction using BEADS with method of residual variances for alpha tuning 
% Input
% y : noisy signal
% alphas : set of relation parameters 
% Output 
% z : denoised signal
% b : baseline

%Compute estimate of redidual variance of noisy signal y by modified one 
% nearest neighbor estimator
x=linspace(1,length(y),length(y))';
[V] = modoneNNestimator(y,x,type);

fc = 0.005;             % fc : cut-off frequency (cycles/sample)
d = 1;                  % d : filter order parameter (d = small positive integer)
Nit = 50;
EPS0 = 1e-5;            % for x itself
EPS1 = 1e-5;            % for derivatives
r = 1;                  % to set the asymmetric ratio
pen = 'L1_v1';          % to choose penalty function

Errs = zeros([1, length(alphas)]); % store residual variances
%Compute residual variances for each hyperparameter value 
for i = 1:length(alphas)
    alpha = alphas(i);
    [f,b,~,~,~] = beads(y, d, fc, r, alpha, pen, Nit, EPS0, EPS1);
    Errs(:,i) = (sum((y-(f+b)).^2))/length(y);
end
% Choose lambda
[~,i] = min(abs(Errs - V));
[z,b,~,~,~] = beads(y, d, fc, r, alphas(i), pen, Nit, EPS0, EPS1);

end