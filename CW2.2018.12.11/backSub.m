function x = backSub(A)
%BACKSUB    Backward Substitution
%   x = backSub(A) performs the backward- 
%   substitution step of Gaussian elimination 
%   without row interchanges to obtain the 
%   solution vector x of A. The input matrix A is 
%   an n by n+1 augmented matrix in echelon form.

% Extract number of rows n from A
n = size(A,1);

% Start Backward Substitution loop
for i = n:-1:1
    
    sum=0; % Initialise the summation of aij*xj
        
    % Obtain summation of aij*xj
    for j = n:-1:i+1
        sum = sum+x(j,1)*A(i,j);
    end
    
    x(i,1) = (A(i,n+1)-sum)/A(i,i);
end