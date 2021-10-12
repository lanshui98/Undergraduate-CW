function [p_vec,fp_vec] = bisect(f,a,b,Nmax)
% This function is designed to find the root of f in interval [a,b] by 
% bisection method with the number of iteraction being Nmacx. 
% It will also produce all Pn and f(Pn) during the iteraction.

for i=1:Nmax
    p=(a+b)/2;
    f(p);
    if f(p)*f(a)>0
        a=p;
    else
        b=p;
    end
    p_vec(i,1)=p;
    fp_vec(i,1)=f(p);
end


