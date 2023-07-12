function [z] = autoCorrSG(y,ws,type)
% denoising using SG with method of autocorrelation for ws tuning 
% Input
% y : noisy signal
% ws : set of window sizes
% Output 
% z : denoised signal

% compute autocorrelation value 
e = diff(y);
if strcmp(type,'mean')
    p = autocorrelation(e,'mean');
elseif strcmp(type,'median')
    p = autocorrelation(e,'percentile');
end

order = 2;

Ps = zeros([1, length(ws)]);
for i = 1:length(ws)
    sgf = sgolayfilt(y,order,ws(i));
    Ps(:,i) = autocorrelation(y-sgf,type);
end
% Choose w
[minval,i] = min(abs((Ps - p)));
z = sgolayfilt(y,order,ws(i));
end