function [A,b] = produce(p,n)

% produce
% [A,b] = produce(p,n) returns a matrix A of size nxn 
% and a vector b of size nx1. The input p is non-positive.

A = (2+p)*diag(ones(n,1));

U = (-1)*diag(ones(n-1,1),1);
L = (-1)*diag(ones(n-1,1),-1);
A = A + U + L;

b = ones(n,1);