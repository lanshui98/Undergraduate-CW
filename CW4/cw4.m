%% G12ISC 2018-2019 Coursework 4
% 
% Subject: Numerical Calculus
% Studen ID: 4336432

clear all
close all
clc


%% Question 1 and Question 2

% Q1:

% Following code produces three approximations to the integral with a
% single interval.

f = @(x) (x.^4)-(x.^3);

% Midpoint rule
MR = 2*(1/4)*f(1/2)

% Trapezoidal rule
TR = (1/2)*(1/2)*(f(1/4)+f(3/4))

% Simpson’s rule
SR = (1/3)*(1/4)*(f(1/4)+4*f(1/2)+f(3/4))

% Q2:

% Following code produces the errors in the three approximations.

% Exact value of the integral
I = @(x) (1/5)*(x.^5)-(1/4)*(x.^4);
I0 = I(3/4)-I(1/4);

% Midpoint rule
error_MR = abs(I0-MR)

% Trapezoidal rule
error_TR = abs(I0-TR)

% Simpson’s rule
error_SR = abs(I0-SR)


%% Question 3

clear all
close all
clc

% Following code produces a table containing approximations to integral 
% using Composite Trapezoidal Rule with different strips and their errors.

% data
a = 0;
T = 128;
f = @(t) 35*(1-exp( -t/(3*T) + (t/(2*T)).^2 - (t/T).^3 ));

% True value
v = 1136.067066175433;

% Construct the vector of m
m = (1:512)';

% Initialize
v_m = zeros(512,1);
error = zeros(512,1);

% Compute the approximations and errors
for i = 1:512
    v_m(i) = CoTrapRule(a,T,f,m(i));
    error(i) = abs(v-v_m(i));
end

% Create table
format long

table(m,v_m,error)

%% Question 4

% Following code produces a plot of the errors using Composite Trapezoidal 
% Rule versus h which is T/m.

% Compute h
h = zeros(512,1); % Initialize
for i = 1:512
    h(i) = T/m(i);
end

% plot
loglog(h,error,'bo')
hold on
plot([1e1,1e3],[1e-3,1e1],'k-')

% Format figure
title('errors versus h using CTR');
xlabel('h');
ylabel('error');
legend('forwElim','O(h^2)','location','bestoutside');

% From the plot, the asymptotic behaviour of error as h tends to 0 is
% O(h^2). And this is verified by the theorem 4.5 which says the error term
% of Composite Trapezoidal Rule is bounded by constant*abs(h^2).


%% Question 5

% Following code produces a rough bound on the errors as a function of m.

% The theorem 4.5 says the absolute value of the error is equal to
% (T-0)/12*(h^2)*f''(u) where u is between 0 and T. Therefore, the error
% bound is equal to (T/12)*(h^2)*max(f''(x)).

syms t
T = 128;
f = 35*(1-exp( -t/(3*T) + (t/(2*T)).^2 - (t/T).^3 ));
df_2 = diff(f,2); % the second derivative of f

% Find the rough maximum value of g
x = (0:T);
max = 0;
for i = 1:T+1
    max = max(vpa(subs(df_2,t,x(i))),max);
end
    
% Bound of error as function of the number of strips m
e = @(m) (T/12)*((T/m)^2)*max;

% When m = 16:
error_m = e(16);
fprintf('The error bound for m = 16 is %f',error_m);

%% Question 7

clear all
close all
clc

% Q3':
% Following code produces a table containing approximations to integral 
% using Composite Midpoint Rule with different strips and their errors.

% data
a = 0;
T = 128;
f = @(t) 35*(1-exp( -t/(3*T) + (t/(2*T)).^2 - (t/T).^3 ));

% True value
v = 1136.067066175433;

% Construct the vector of m
m = (1:512)';

% Initialize
v_m = zeros(512,1);
error = zeros(512,1);

% Compute the approximations and errors
for i = 1:512
    v_m(i) = CoMidRule(a,T,f,m(i));
    error(i) = abs(v-v_m(i));
end

% Create table
format long

table(m,v_m,error)

% Q4':
% Following code produces a plot of the errors using Composite Midpoint 
% Rule versus h which is T/(2*m).

% Compute h
h = zeros(512,1); % Initialize
for i = 1:512
    h(i) = T/(2*m(i));
end

% plot
loglog(h,error,'bo')
hold on
plot([1e0,1e2],[1e-1,1e3],'k-')

% Format figure
title('errors versus h using CMR');
xlabel('h');
ylabel('error');
legend('forwElim','O(h^2)','location','bestoutside');

% From the plot, the asymptotic behaviour of error as h tends to 0 is
% O(h^2). And this is verified by the lemma of theorem 4.5 which says 
% the error term of Composite Midpoint Rule is bounded by constant*abs(h^2).

%% Question 9

clear all
close all
clc

% Q3'':
% Following code produces a table containing approximations to integral 
% using Composite Simpson's Rule with different strips and their errors.

% data
a = 0;
T = 128;
f = @(t) 35*(1-exp( -t/(3*T) + (t/(2*T)).^2 - (t/T).^3 ));

% True value
v = 1136.067066175433;

% Construct the vector of m
m = (1:512)';

% Initialize
v_m = zeros(512,1);
error = zeros(512,1);

% Compute the approximations and errors
for i = 1:512
    v_m(i) = CoSimRule(a,T,f,m(i));
    error(i) = abs(v-v_m(i));
end

% Create table
format long

table(m,v_m,error)

% Q4'':
% Following code produces a plot of the errors using Composite Simpson's 
% Rule versus h which is T/(2*m).

% Compute h
h = zeros(512,1); % Initialize
for i = 1:512
    h(i) = T/(2*m(i));
end

% plot
loglog(h,error,'bo')
hold on
plot([1e-1,1e1],[1e-5,1e3],'k-')

% Format figure
title('errors versus h using CSR');
xlabel('h');
ylabel('error');
legend('forwElim','O(h^4)','location','bestoutside');

% From the plot, the asymptotic behaviour of error as h tends to 0 is
% O(h^4). And this is verified by the lemma of theorem 4.5 which says 
% the error term of Composite Simpson's Rule is bounded by constant*abs(h^4).

%% Question 11

clear all
close all
clc

% Following code produces the plot of v(t) versus t = 0,1,2...128. And the 
% error in approximate v(t) for each t is less than 0.05.

t = (0:218); % Initialize t

% data
a = 0;
T = 128;

% The second derivative of f
syms t
f = 35*(1-exp( -t/(3*T) + (t/(2*T)).^2 - (t/T).^3 ));
df2 = diff(f,2); 

% Compute df2_max for each t
df2_max = zeros(219,1); % Initialize df2_max
for j = 1:T+1
    x = (0:t(j));
    df2_max(j) = 0;
    for i = 1:t(j)+1
        df2_max(j) = max(df2_max(j),(vpa(subs(df2(j),t,x(i)))));
    end
end

% Compute the approximate v(t) for each t such that the error is less than
% 0.05
v = zeros(219,1); % Initialize v(t)
for i = 1:219
    TOL = 0.05;
    f = @(t) 35*(1-exp( -t/(3*T) + (t/(2*T)).^2 - (t/T).^3 ));
    v(i) = VelTrapRule(i-1,TOL,f,df2_max(i));
end

% plot
plot(v,t);
% Format figure
title('v(t) versus t = 0,1,2,...,218');
xlabel('t');
ylabel('v(t)');

% Following code produces the plot of m versus the t where m is the 
% appropriate number of strips used for approximation for each t.

% Compute the m for each t
e = @(m) (t/12)*((t/m)^2)*df2_max;
m = zeros(219,1); % Initialize m
for j = 1:219
    for i = 1:512
        if e(i)<TOL
            m(j) = i;
            break
        end
    end
end

% plot
plot(m,t); 
% Format figure
title('appropriate number of strips versus t = 0,1,2,...,218');
xlabel('t');
ylabel('m');

