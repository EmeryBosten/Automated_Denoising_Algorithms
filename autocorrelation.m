function [p] = autocorrelation(e,type)
% autocorrelation function to compute correlation of vector
% noise vector (can be residual vector)
% p = 0 white noise
% p > 0 red noise
% p < 0 blue noise
% different types of computations to render more robust by using 
% median, L1 norm or percentile of data or combination thereof
% based on 'Automatic Selection of Optimal Savitzky-Golay Smoothing' by
% Gabriel Vivo´-Truyols and Peter J. Schoenmakers
if strcmp(type,'median')
    p = 1 - (median(diff(e).^2)/(2*median(e.^2)));
elseif strcmp(type,'L1')
    p = 1 - (sum(abs(diff(e)))/(2*sum(abs(e))));
elseif strcmp(type,'percentile')
    e = e(e<prctile(e,90) & e > prctile(e,10));
    e(e==0) = [];
    p = 1 - (sum(diff(e).^2)/(2*sum(e.^2)));
elseif strcmp(type,'percmed')
    e = e(e<prctile(e,90) & e > prctile(e,10));
    e(e==0) = [];
    p = 1 - (median(diff(e).^2)/(2*median(e.^2)));
elseif  strcmp(type,'mean')
    p = 1 - (sum(diff(e).^2)/(2*sum(e.^2)));
end