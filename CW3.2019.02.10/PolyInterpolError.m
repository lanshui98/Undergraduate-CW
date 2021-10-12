function e = PolyInterpolError(a,b,f,z) 
% PolyInterpolError   Polynomial interpolant error.
%     e = PolyInterpolError(a,b,f,z)  computes the max-norm error e
%           of the n-th degree polynomial interpolating f 
%           at n+1 equally-spaced points in the interval [a,b] contained
%           in vector z

% Define (suffient number of) evaluation points in [a,b]
x = GridEq(160,a,b);
    
% Compute the max-norm error in [a,b] between the interpolant and f

e = 0; % Initialize error

for i = 1:160
    e = max(abs(f(x(i))-PolynInterp(f,x(i),z)),e);
end
