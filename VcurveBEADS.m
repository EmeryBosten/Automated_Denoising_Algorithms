function [z,b] = VcurveBEADS(y,alphas)
% denoising and baseline correction using BEADS with L- and V- curve for alpha tuning 
% Input
% y : noisy signal
% alphas : set of relation parameters (determines relation between regularization parameters lambda)
% Output 
% z : denoised signal
% b : baseline
% hyperparameters to be discussed - see BEADS paper 
fc = 0.006;             % fc : cut-off frequency (cycles/sample)
d = 1;                  % d : filter order parameter (d = small positive integer)

% defined hyperparameters
Nit = 50;
EPS0 = 1e-5;            % for x itself
EPS1 = 1e-5;            % for derivatives
r = 1;                  % to set the asymmetric ratio
pen = 'L1_v1';          % to choose penalty function

Ps = zeros([1 length(alphas)]); % store penalties for each lam
Rs = zeros([1 length(alphas)]); % store residuals for each lam
xs = zeros([length(y) length(alphas)]); % store estimated for each lam 
fs = zeros([length(y) length(alphas)]); % store baseline for each lam

% L-curve
for i = 1:length(alphas)
    alpha = alphas(i);
    [x,f,~,R,P] = beads(y, d, fc, r, alpha, pen, Nit, EPS0, EPS1);
    Rs(i) = R;
    Ps(i) = P;
    xs(:,i) = x;
    fs(:,i) = f;
end

% V-Curve
% Implementation follows "L- and V-curves for optimal smoothing"
%by Gianluca Frasso and Paul HC Eilers
 d3Mat = [log10(Rs) ; log10(Ps) ; alphas]';
 d3MatB = sortrows(d3Mat,1);

% Compute V curve 
distance = zeros([1 length(alphas)-1]);
alpha = zeros([1 length(alphas)-1]);
for i = 1:length(alphas)-1
    distance(i) = sqrt((d3MatB(i+1,2)-d3MatB(i,2))^2 + (d3MatB(i+1,1)-d3MatB(i,1))^2);
    alpha(i) = (d3MatB(i+1,3) + d3MatB(i,3))/2;
end 

%find minimum distance
[M1,I] = min(distance);
% run BEADS with optimal lambda
[z,b,~,~,~] = beads(y, d, fc, r, alpha(I), pen, Nit, EPS0, EPS1);
end