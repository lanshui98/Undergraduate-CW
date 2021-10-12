function [B,P] = forwElimPPperm(A)

%forwElimPPperm
%   [B,P] = forwElimPPperm(A) performs the forward- 
%   elimination step of Gaussian elimination 
%   with row interchanges to obtain the echelon 
%   form B of A by partial pivoting as well as the 
%   permutation matrix P corresponding to the pivoting process.
%   The input matrix A is an n by n+1 augmented matrix.

% Extract number of rows n from A
n = size(A,1);

% Initialize the max and index
max=0;
index=0;

% LU Factorization when Gaussian Elimination works with Partial Pivoting 
% has the form PA = LU where P is the permutation matrix. In order to
% obtain P, we can first calculate L and U.
% In addition, U equals the echelon form B of A.

% Initialize the permutation matrix P which is an identity matrix.
P=eye(n);

% Start Gaussian elimination loop
for i = 1:n-1
    
    % Find the maximum value in the column under the current checked element
    % and return its row position
    for l = i:n
        if abs(A(l,i))>max
            max = abs(A(l,i));
            index = l;
        end
    end
    
    % Perform the row exchanges in both A and P which are the same.
    A([i index],:)=A([index i],:);
    P([i index],:)=P([index i],:);
    
    % Eliminate column i
    for j = i+1:n
        
        % Compute multiplier
        m = A(j,i)/A(i,i);
        
        % Replace Ej by Ej-m*Ei
        A(j,i) = 0;
        for k = i+1:n+1
          A(j,k) = A(j,k) - m*A(i,k);
        end
        
    end
    
end

% Return one of the outputs B
B = A;







