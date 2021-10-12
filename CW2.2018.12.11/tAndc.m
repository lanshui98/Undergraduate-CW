function [T,c] = tAndc(A,b)

% tAndc
% [T,c] = tAndc(A,b) returns T and c of 
% linear system Ax = b under the Jacobi method.

D = diag(diag(A));
L = -1*tril(A,-1);
U = -1*triu(A,1);

% Outputs
T = D\(L+U);
c = D\b;