function [x_mat,Niter,resNorm_vec] = genIterTol(T,c,x0,Nmax,tol,A,b)

%genIterTol

% [x_mat,Niter,resNorm_vec] = genIterTol(T,c,x0,Nmax,tol,A,b) performs
% the general iteration method as genIter.m but with a stopping criterion 
% based on the l2 norm of the residual vector. 
% The input consists additionally of a real number tol, matrix A and vector b.
% The output consists additionally of the total number of iterations Niter
% required to reach the stopping criterion (or to reach Nmax, whichever
% comes first) and resNorm_vec which stores all residual vectors.

% Put in the initial approximation vector
x_mat(:,1) = x0;

% Initialize the Niter
Niter = 0;

% Start the iteration x(k)=T*x(k−1)+c
for i=2:Nmax+1
    x0 = T*x0+c;
    error = b-A*x0; % the residual vector
    x_mat(:,i) = x0;
    Niter = Niter+1;
    resNorm_vec(i-1,1) = norm(error,2);
    if norm(error,2) < tol
        break
    end
end

