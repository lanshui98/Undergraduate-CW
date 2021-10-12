function x_mat = genIter(T,c,x0,Nmax)

%genIter
% x_mat = genIter(T,c,x0,Nmax) performs the iteration method x(k)=T*x(k−1)+c
% to obtain a matrix x_mat that stores all approximation vectors x(k). 
% The input consists of an iteration matrix T, an vector c, the initial 
% approximation vector x0 and the total number of iterations to be peformed Nmax.

% Extract number of rows n of c
n=length(c);

% Initialize x_mat
x_mat = zeros(n,Nmax+1);

% Put in the initial approximation vector
x_mat(:,1) = x0;

% Start the iteration x(k)=T*x(k−1)+c
for i=2:Nmax+1
    x0 = T*x0+c;
    x_mat(:,i) = x0;
end
    
