function [z] = resVarSG(y,ws,type)
% denoising using SG smoother with method of residual variances for window size tuning 
% Input
% y : noisy signal
% ws : set of window sizes 
% Output 
% z : denoised signal

%Compute estimate of redidual variance of noisy signal y by modified one 
% nearest neighbor estimator
x=linspace(1,length(y),length(y))';
[V] = modoneNNestimator(y,x,type);
order = 2;

Errs = zeros([1, length(ws)]);
for i = 1:length(ws)
    sgf = sgolayfilt(y,order,ws(i));
    Errs(:,i) = (sum((y-sgf).^2))/length(y);
end
% Choose w
[~,i] = min(abs(Errs - V));
z = sgolayfilt(y,order,ws(i));
end