function B = forwElimPP(A)

%FORWELIMPP    Forward Elimination with Partial Pivoting
%   B = forwElimPP(A) performs the forward- 
%   elimination step of Gaussian elimination 
%   with row interchanges to obtain the 
%   echelon form B of A by partial pivoting. 
%   The input matrix A is an n by n+1 augmented matrix.

% Extract number of rows n from A
n = size(A,1);

% Initialize the max and index
max=0;
index=0;

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
    
    %perform row exchange in the nrow vector
    A([i index],:)=A([index i],:);
    
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

% Return output B
B = A;