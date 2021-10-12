%% G12ISC 2018-2019 Coursework 5
% Student ID: 4336432
% Subject: Numerical ODE's


clear all
close all
clc
 
%% Question 1

% Following code produces a table and a figure comparing the approximations
% by Euler's method with exact solutions.

f = @(t,y) y-t^2+1; % y'
y = @(t) (t+1)^2-0.5*exp(t); % exact solution

% Data
a = 0;
b = 3;
h = 0.2;
w0 = 0.5;
N = (b-a)/h;

% Create table
format long g
[t,w]=Euler(a,b,N,w0,f);
Y = zeros(N+1,1);
error = zeros(N+1,1);
for i = 1:N+1
    Y(i) = y(t(i));
    error(i) = abs(w(i)-Y(i));
end
table(t,w,Y,error)

% Create figure
plot(t,w)
hold on
plot(t,Y)
% Format figure
title('Approximations by Euler method and exact solutions of IVP')
legend('approximation wi','exact solution Yi','location','bestoutside')
xlabel('ti')

%% Question 2

clear all
close all
clc

% Following code produces a table and a figure to invesgate the asymptotic 
% rate of convergence of Euler's method.

f = @(t,y) y-t^2+1; % y'
y = @(t) (t+1)^2-0.5*exp(t); % exact solution

% Data
a = 0;
b = 3;
w0 = 0.5;

h = zeros(7,1); % initialise
e = zeros(7,1);

% Create table
format long g
for i = 2:8
    h(i-1) = 1/(2^i);
    N = (b-a)/h(i-1);
    [t,w] = Euler(a,b,N,w0,f);
    e(i-1) = abs(w(N+1)-y(3));
end
table(h,e)

% Create figure
loglog(h,e)
% Format figure
title('error of Euler Method versus h')
xlabel('h')
ylabel('error')
grid on

% Therefore, from the plot, errors of Euler's Method converge linearly.

%% Question 4

clear all
close all
clc

% Following code produces a table and a figure comparing the approximations
% by Modified Euler's method with exact solutions.

f = @(t,y) y-t^2+1; % y'
y = @(t) (t+1)^2-0.5*exp(t); % exact solution

% Data
a = 0;
b = 3;
h = 0.2;
w0 = 0.5;
N = (b-a)/h;

% Create table
[t,w]=ModifiedEuler(a,b,N,w0,f);
Y = zeros(N+1,1);
error = zeros(N+1,1);
for i = 1:N+1
    Y(i) = y(t(i));
    error(i) = abs(w(i)-Y(i));
end
table(t,w,Y,error)

% Create figure
format long g
plot(t,w)
hold on
plot(t,Y)
% Format figure
title('Approximations Modified Euler method and exact solutions of IVP')
legend('approximation wi','exact solution Yi','location','bestoutside')
xlabel('ti')

%% Question 5

clear all
close all
clc

% Following code produces a table and a figure to invesgate the asymptotic 
% rate of convergence of Modified Euler's method.

f = @(t,y) y-t^2+1; % y'
y = @(t) (t+1)^2-0.5*exp(t); % exact solution

% Data
a = 0;
b = 3;
w0 = 0.5;

h = zeros(7,1); % initialise
e = zeros(7,1);

% Create table
format long g
for i = 2:8
    h(i-1) = 1/(2^i);
    N = (b-a)/h(i-1);
    [t,w] = ModifiedEuler(a,b,N,w0,f);
    e(i-1) = abs(w(N+1)-y(3));
end
table(h,e)

% Create figure
loglog(h,e)
% Format figure
title('error of Modified Euler Method versus h')
xlabel('h')
ylabel('error')
grid on

% Therefore, from the plot, errors of Modified Euler's Method converge 
% quadratically.

%% Question 7

clear all
close all
clc

% Following code produces a table and a figure comparing the approximations
% by Runge–Kutta Order Four method with exact solutions.

f = @(t,y) y-t^2+1; % y'
y = @(t) (t+1)^2-0.5*exp(t); % exact solution

% Data
a = 0;
b = 3;
h = 0.2;
w0 = 0.5;
N = (b-a)/h;

% Create table
format long g
[t,w]=RuKuMeth(a,b,N,w0,f);
Y = zeros(N+1,1);
error = zeros(N+1,1);
for i = 1:N+1
    Y(i) = y(t(i));
    error(i) = abs(w(i)-Y(i));
end
table(t,w,Y,error)

% Create figure
plot(t,w)
hold on
plot(t,Y)
% Format figure
title('Approximations by Runge–Kutta Order Four method and exact solutions of IVP')
legend('approximation wi','exact solution Yi','location','bestoutside')
xlabel('ti')


%% Question 8

clear all
close all
clc

% Following code produces a table and a figure to invesgate the asymptotic 
% rate of convergence of Runge–Kutta Order Four method.

f = @(t,y) y-t^2+1; % y'
y = @(t) (t+1)^2-0.5*exp(t); % exact solution

% Data
a = 0;
b = 3;
w0 = 0.5;

h = zeros(7,1); % initialise
e = zeros(7,1);

% Create table
format long g
for i = 2:8
    h(i-1) = 1/(2^i);
    N = (b-a)/h(i-1);
    [t,w] = RuKuMeth(a,b,N,w0,f);
    e(i-1) = abs(w(N+1)-y(3));
end
table(h,e)

% Create figure
loglog(h,e)
% Format figure
title('error of Runge–Kutta Order Four method versus h')
xlabel('h')
ylabel('error')
grid on

% Therefore, from the plot, the order of convergence of this method is
% four.

%% Question 10

clear all
close all
clc

% Following code produces a plot of Euler's approximations w1,i and the
% exact solutions and figure out for what value of N will the computed
% approximation become completely unreliable.

T = 7;
N = 3000;
theta = 2^(1/2);
alpha = pi/10;
delta = 0;

[t,w1,w2]=EulerSys(T,N,theta,alpha,delta);
y = @(t) alpha*cos(theta*t); % exact solution
Y = zeros(N+1,1); % initialise

% Create figure
for i=1:N+1
    Y(i) = y(t(i));
end
plot(t,w1)
hold on
plot(t,Y)
% Format figure
title('approximations w1,i and the exact solutions')
xlabel('ti')
legend('approximation w1,i','exact solution Yi','location','bestoutside')

%% Question 10 cont.

clear all
close all
clc

% That computed approximation become completely unreliable mathematically
% means the approximation go far away from the exact solution. Here 10
% values of n are tested by ploting approximations under each n.

T = 7;
theta = 2^(1/2);
alpha = pi/10;
delta = 0;

% Create figure
n = linspace(10,100,10);
for i = 1:10
    [ti,w1,w2]=EulerSys(T,n(i),theta,alpha,delta);
    plot(ti,w1)
    hold on
end
t = linspace(0,T,100);
y = alpha*cos(theta*t); % exact solution
plot(t,y,'LineWidth',2)
hold off
% Format figure
title('approximations under different n')
xlabel('ti')
legend('n=10','n=20','n=30','n=40','n=50','n=60','n=70','n=80','n=90','n=100','exact solutions','location','bestoutside')

% From the plot, we can see as n becomes smaller, the approximations are
% less accurate. When n is less than 20, the tendency of the approximation
% goes against the exact solution which is a periodic function. Therefore,
% roughly speaking, when n is less than 20, the computed approximation will 
% become completely unreliable.




