function [A, B] = BAfilt(d, fc, N)
% [A, B] = ABfilt(d, fc, N)
%
% Banded matrices for zero-phase high-pass filter.
% The matrices are 'sparse' data type in MATLAB.
%
% INPUT
%   d  : degree of filter is 2d
%   fc : cut-off frequency (normalized frequency, 0 < fc < 0.5)
%   N  : length of signal


b1 = [1 -1];
for i = 1:d-1
    b1 = conv(b1, [-1 2 -1]);
end
b = conv(b1, [-1 1]);

omc = 2*pi*fc;
t = ((1-cos(omc))/(1+cos(omc)))^d;

a = 1;
for i = 1:d
    a = conv(a,[1 2 1]);
end
a = b + t*a;
A = spdiags( a(ones(N,1), :) , -d:d, N, N);     % A: Symmetric banded matrix

B = spdiags(b(ones(N,1), :) , -d:d, N, N);       % B1: banded matrix

e = ones(N, 1);

