function I=CoSimRule(a,b,f,m)

% CoSimRule   Composite Simpson's Rule
%   I=CSR(a,b,f,m) approximates the integral of f(x) between 
%   a and b using the Composite Simpson's Rule with m strips (i.e. via 2m+1 points)
%     f is a function handle
%     a and b are scalars
%     m is an integer greater or equal than 1


% Define h
h=(b-a)/(2*m); 

% Define points to evaluate function
x=zeros(2*m+1,1);  
for i=0:2*m
    x(i+1)=a+i*h;
end

% Compute the approximations via composite simpson's rule:
I=0;
for i=1:m
    c = 2*i-1;
    I=I+f(x(c))+4*f(x(c+1))+f(x(c+2));
end
I=h/3*I;



