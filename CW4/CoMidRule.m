function I=CoMidRule(a,b,f,m)

% CoMidRule   Composite Midpoint Rule
%   I=CMR(a,b,f,m) approximates the integral of f(x) between 
%   a and b using the Composite Midpoint Rule with m strips (i.e. via m points)
%     f is a function handle
%     a and b are scalars
%     m is an integer greater or equal than 1


% Define h
h=(b-a)/(2*m); 

% Define points to evaluate function
x=linspace(a+h,b-h,m);

% Compute the approximations via composite midpoint rule:
I = 2*h*sum(f(x(1:end)));
