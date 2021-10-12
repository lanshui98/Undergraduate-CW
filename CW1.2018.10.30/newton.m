function [p_vec] = newton(f,df,p0,Nmax,tol)
% This function is designed to find the root of f by Newton's method with df being f',
% p0 being the initial guess, Nmacx being the maximum number of iteraction and tol being the tolerance. 
% It will produce all Pn during the iteraction.

for i=1:Nmax
    p=p0-f(p0)/df(p0);
    p_vec(i,1)=p;
    if abs(p-p0)<tol
        break;
    else
        p0=p;
    end
end
