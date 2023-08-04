function z = WSsparse(y, lambda,d)
% Whittaker smoother with use of sparsity
% Implementation follows "A perfect smoother" by P. Eilers
m = length(y);
E = speye(m);
D = diff(E,d);
C = chol(E + lambda * (D' * D));
z = C\(C'\y');
end