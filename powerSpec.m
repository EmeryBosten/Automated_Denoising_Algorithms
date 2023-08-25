function [power_spectrum, f] = powerSpec(signal,fs)
%function that computes power spectrum of signal
% INPUT 
% signal = signal
% fs = sampling frequency
% OUTPUT
% power_spectrum = corresponding power spectrum of signal
% f = frequency values 
N = length(signal);
fft_result = fft(signal);

% Calculate the corresponding frequency values
f = (0:N-1)*(fs/N);

% Calculate the power spectrum (magnitude squared of FFT result)
power_spectrum = abs(fft_result).^2;
end