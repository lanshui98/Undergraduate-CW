function I=CoTrapRule(a,b,f,m)
% CoTrapRule   Composite Trapezoidal Rule
%   I = CTR(a, b,f, m) approximates the integral of f(x) between 
%   a and b using the Composite Trapezoidal Rule with m strips (i.e. via m+1 points)
%     f is a function handle
%     a and b are scalars
%     m is an integer greater or equal than 1

h=(b-a)/m; %define strip length / spacing between points


%%%% version with loops

%%%Define points to evaluate function (at end points of each strip) 
x=zeros(m+1,1);  
for i=0:m
    x(i+1)=a+i*h;
end

%%%Compute the approximations via composite trapezoidal rule:
I=0;
for i=2:m
    I=I+f(x(i));
end
I=h/2*(2*I+f(x(1))+f(x(m+1)));


%%%Alternative shorter version of the code:
%h=(b-a)/m;
%x=linspace(a,b,m+1);
%I=h/2*(f(x(1))+f(x(end))+2*sum(f(x(2:end-1))));
