function y = PolynInterp(f,x,z)
%PolyInterp   Polynomial interpolant.
%     y = PolynInterp(f,x,z) computes the interpolating
%     polynomial P such that P(z) = f(z), and evaluates P in 
%     the points given by x, thereby returning y=P(x).


% Get the number of points from the input z
n = size(z,2)-1;

% Initialize the output y 
m = size(x,2);
y = zeros(1,m);

% Compute the output y using a loop

for i = 1:m
    for j = 0:n
        y(i) = y(i)+f(z(j+1))*LagrPolyn(j,x(i),z);
    end
end    
