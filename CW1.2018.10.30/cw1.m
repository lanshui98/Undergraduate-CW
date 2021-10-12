%% MATH2019 / G12ISC (2018-2019) Coursework 1
% NAME: LAN SHUI   STUDENT ID: 4336432

clear all;
close all;
clc

%% Question 2

% The code produces a table containing n, pn and f(pn) during the bisection method.

format long
f=@(x)x.^3+4*x.^2-10;
[p_vec,fp_vec]=bisect(f,1,2,20);
n=[1:20]';
table(n,p_vec,fp_vec)


%% Question 3

% The code creates a figure ploting the error versus n during the bisection method.

p=1.365230013;
for i=1:20
    e=abs(p-p_vec(i));
    semilogy(i,e,'bo')
    hold on
end

%% Question 5

% The code produces a table containing n and pn during the fixed- point iteration.

clear all;
close all;
clc

format long
f=@(x)x.^3+4*x.^2-10;
g=@(x)x-1/12*f(x);
p_vec=fpiter(g,1,20);
n=[1:20]';
table(n,p_vec)

%% Question 6

% The code creates a figure ploting the error versus n during the fixed- point iteration.

p=1.365230013;
for i=1:20
    e=abs(p-p_vec(i));
    semilogy(i,e,'bo')
    hold on
end

%% Question 7

% This code produces the pmax that for any p0 >= pmax, the fixed-point iteration won't converge.

% A general fixed-point method converges linearly if it converges. So a
% criterion that can be taken to judge if the iteration converges is the 
% difference between p and p20 after 20 iterations. If it's large, say large
% than 0.001, the iteration doesn't converge.

% This experimental method isn't accurate enough as 0.001 and 20 don't
% have a very strong theoretical basis. The accuracy is also influenced by
% the amount of data being chosen.

a=1.37; 
% The experimental data {an} being chosen is an arithmetic progression with
% a0=1.37, d=0.5 and n=10.
for i=1:10
    p_vec=fpiter(g,a,20);
    if abs(p-p_vec(20,1))>0.001
        pmax=a
        break;
    else
        a=a+0.5;
    end
end

%% Question 9

% This code produces a table containing n and pn during Newton's method.

clear all;
close all;
clc

format long
f=@(x)x.^3+4*x.^2-10;
df=@(x)3*x.^2+8*x;
p_vec=newton(f,df,1,40,1e-10);
n=[1:size(p_vec,1)]';
table(n,p_vec)

% The code creates a figure ploting the error versus n during Newton's method.

p=1.365230013;
for i=1:size(p_vec,1)
    e=abs(p-p_vec(i));
    semilogy(i,e,'bo')
    hold on
end

%% Question 10

% This part invesgates the convergence behaviour of Newton's method for values
% of p0 that approach 0 from above.

a=0.1; 
% The experimental data {an} being chosen is an arithmetic progression with
% a0=0.1, d=0.5 and n=10.

% The code produces figures ploting pn versus n during the iteration for each p0.

for i=1:10
    p_vec=newton(f,df,a,40,1e-10);
    n=size(p_vec,1);
    subplot(2,5,i);
    for j=1:n
        semilogy(j,p_vec(j),'bo')
        hold on 
    end
    a=a+0.5;
end
        
% From the gragh, it can be figured out that the Newton's method converges
% quite quickly for p0 near 0 from above.
