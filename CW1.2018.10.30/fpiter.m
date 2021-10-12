function [p_vec] = fpiter(g,p0,Nmax)
% This function is designed to find the root of f by fixed- point iteration 
% with g being a function handle, p0 being the initial guess and Nmacx being the number of
% iteraction. 
% It will also produce all pn during the iteraction in a column vector.

for i=1:Nmax
    p=g(p0);
    p0=p;
    p_vec(i,1)=p;
end