function [A,b] = testMat(n)
%TESTMAT   Test matrix and vector
%   [A,b] = testMat(n) returns 
%   a matrix A of size nxn 
%   and a vector b of size nx1

A = diag(4*ones(n,1));
for i=1:n-1
    U = diag(1/2^i*ones(n-i,1),i);
    L = diag(1/2^i*ones(n-i,1),-i);
    A = A + U + L;
end

b = ones(n,1);

