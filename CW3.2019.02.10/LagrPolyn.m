function y = LagrPolyn(k,x,z)

% LagrPoly   Lagrange polynomial.
%     y = LagrPoly(k,x,z) computes the nth Lagrange polynomial 
%     based on the (n+1) points in the vector z, and evaluates the 
%     Lagrange polynomial in the points given by the vector x, 
%     thereby returning the vector y=L_k(x).


% Get the number of points from the input z
n = size(z,2)-1;

% Initialize the output y 
m = size(x,2);
y = ones(1,m);
    
% Compute the output y using a loop

for i = 1:m
    for j = 0:n
        if j ~= k
            y(i) = y(i)*(x(i)-z(j+1))/(z(k+1)-z(j+1));
        end
    end
end

