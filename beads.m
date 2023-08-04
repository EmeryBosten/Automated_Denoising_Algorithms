function [x,f,cost,R,P] = beads(y, d, fc, r, alpha, pen, Nit, EPS0, EPS1)
% BEADS algorithm following "Chromatogram baseline estimation and denoising using
% sparsity (BEADS)" by Xiaoran Ning, Ivan Selesnick, and Laurent Duval
% [x,f,cost] = beads(y, d, fc, r, lam0, lam1, lam2, pen, Nit, EPS0, EPS1)
% Chromatogram baseline estimation and denoising using sparsity (BEADS)
% INPUT
%       y: Noisy observation
%       d: Order of the high-pass filter H
%       fc: Normalized cut-off frequency of the high-pass filter H
%       r: Aymmetric ratio
%       lam0, lam1, lam2: Regularization parameters -> has been summarized 
%       by a single relation parameter alpha
%       pen  : 'L1_v1' or 'L1_v2'
%       Nit: Number of iterations
%       EPS0, EPS1: small numbers
%       BA_filt: A function used to generate matrices A and B such that H = B*A^{-1}
% OUTPUT
%       x: Estimated sparse derivative signal component
%       f: Estimated low-pass component
%       cost: Cost function history
%       R : Residuals (Fidelity term)
%       P : Penalty (Penalty term)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch pen
    case 'L1_v1'
        phi = @(x) sqrt(abs(x).^2 + EPS1);
        wfun = @(x) 1./(sqrt(abs(x).^2 + EPS1));
    case 'L1_v2'
        phi = @(x) abs(x) - EPS1 * log(abs(x) + EPS1);
        wfun = @(x) 1./( abs(x) + EPS1);
        %     case 'log'
        %         phi = @(x) (1/a)*log(1 + a * sqrt(abs(x).^2 + EPS1));
        %         wfun = @(x) 1./ ( sqrt(abs(x).^2 + EPS1) + a * (abs(x).^2 + EPS1));
    otherwise
        disp('penalty must be L1_v1, L1_v2, or log')
        x = []; cost = []; f = [];
        return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lam0 = 0.1*alpha;
lam1 = 5*alpha;
lam2 = 5*alpha;

y = y(:);
x = y;
cost = zeros(1,Nit);
N = length(y);
[A, B] = BAfilt(d, fc, N);
H = @(x) B*(A\x);
e = ones(N-1,1);
D1 = spdiags([-e e],[0 1],N-1,N);
D2 = spdiags([e -2*e e], 0:2, N-2, N);
D = [D1;  D2];
w = [lam1 * ones(N-1,1); lam2 * ones(N-2,1)];


b = (1-r)/2 * ones(N,1);

d = B'*B * (A\y) - lam0 * A'*b;


for i = 1:Nit
    
    
    
    temp = ((1 + r)/4)./ (abs(EPS0) * ones(N,1));
    temp(abs(x)>EPS0) = ((1 + r)/4) ./  abs(x(abs(x)>EPS0));
    
    Lambda = spdiags( w.*wfun(D*x), 0, 2*N-3, 2*N-3);
    
    M = 2 * lam0 * spdiags( temp, 0, N, N) + D'*Lambda*D;
    Q = B'*B + A'*M*A;
    
    x = A * (Q\d);
    
    R = 0.5 * sum(abs(H(y - x)).^2);
    P =  ( sum(x(x>EPS0)) + (-r) * sum(x(x<-EPS0)) +...
        sum( (1+r)/4/abs(EPS0)*x(abs(x)<=EPS0).^2 + (1-r)/2 * x(abs(x)<=EPS0) + abs(EPS0)*(1+r)/4) )...
        + sum(phi(diff(x)))...
        + sum(phi(diff(x,2)));
    
    cost(i) = 0.5 * sum(abs(H(y - x)).^2)+ ...
        lam0 * ( sum(x(x>EPS0)) + (-r) * sum(x(x<-EPS0)) +...
        sum( (1+r)/4/abs(EPS0)*x(abs(x)<=EPS0).^2 + (1-r)/2 * x(abs(x)<=EPS0) + abs(EPS0)*(1+r)/4) )...
        + lam1 * sum(phi(diff(x)))...
        + lam2 * sum(phi(diff(x,2)));
end

f = y - x - H(y-x);






