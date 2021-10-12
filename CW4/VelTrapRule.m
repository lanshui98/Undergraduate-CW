function vt = VelTrapRule(t,TOL,f,df2_max)

% VelTrapRule
%   vt = VelTrapRule(t,TOL,f,df2_max) approximates the integral of f between 
%   0 and t using the Composite Trapezoidal Rule such that the error of the
%   integral is less than the given tolerance which is TOL.
%     t is a value between 0 and T
%     TOL is a scalar
%     f is a function handle
%     df2_max is a bound for abs(f''(t))

% Compute the bound on the error obtained with CoTrapRule.m as a function
% of the number of strips m:
e = @(m) (t/12)*((t/m)^2)*df2_max;

% Obtain the appropriate value of m such that the error of is less than TOL
for i = 1:512
    if e(i)<TOL
        m = i;
        break
    end
end

% Compute the approximation with obtained m strips
vt = CoTrapRule(0,t,f,m);

