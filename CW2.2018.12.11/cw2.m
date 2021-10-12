%% MATH2019 / G12ISC (2018-2019) Coursework 2
% NAME: Shui Lan,  STUDENT ID: 4336432

clear all;
close all;
clc



%% Question 2

% The code computes the echelon
% form and the solution in case 
% the system of linear equations 
% has an augmented matrix A1.


A1 = [...
     1    -1     2    -1    -8
     1     1     1     0    -2
     2    -2     3    -3   -20
     1    -1     4     3     4
    ];
B1 = forwElim(A1)
x = backSub(B1)

%% Question 3

% The code computes the time forwElim.m takes to obtain 
% the echelon form of a n*n+1 augmented matrix and the time 
% backSub.m takes to obtain the solution vector.

% It also plots the time they take versus n.

clear all;
close all;
clc

n=50;
for i=1:20

    % Obtain the augmented matrix
    [A,b] = testMat(n);
    B = [A,b];

    % compute the time (t1) forwElim.m takes
    tic;

    C = forwElim (B);

    t1(i)=toc;

    % compute the time (t2) backSub.m takes
    tic;

    x = backSub(C);

    t2(i)=toc;
    
    n=n+50;

end

% Plot the computed time for forwElim.m versus n and
% computed time for backSub.m in the same figure

x=linspace(50,1000,20);
loglog(x,t1,'ro');

hold on

loglog(x,t2,'bo');
title('Time taken to compute forwElim.m and backSub.m');
xlabel('n');
ylabel('time t');

legend('forwElim.m','backSub.m');

% What observed from the plots is that the time taken by forwElim.m
% increases as the size of augmented matrix increases, but the time 
% taken by backSub.m dosen't change much as the size of augmented matrix 
% increases.

%% Question 5

% The code computes the echelon
% form with row interchanges by 
% partial pivoting and the solution 
% in case the system of linear equations 
% has an augmented matrix A2.

clear all;
close all;
clc

A2 = [...
     0    -1     2    -1    -1
     1    -1     1     0    -8
     2    -2     3    -3   -20
     1    -1     4     3     4 
    ];
B2 = forwElimPP(A2)
x = backSub(B2)


%% Question 7

clear all;
close all;
clc

% The code computes the permutation matrix P of A3 and the resulting 
% echelon form B by applying forwElim.m to A (= P*A3).

A3 = [...
    26   -34    31   -48    25   -18    13   -32   -28    13     29
     0     0    25    42   -27    31   -14    23   -13   -48     67
    24    39   -38    15   -44    29     0   -13   -41    41   -152
   -39     2     3    43    27     0   -28    34    14     0    362
    18    20   -17   -34    17     1    15    23   -32    25    292
    -4   -35     5    42    22    14    10     7   -45    31    281
     0    45   -10    29    14    45   -11   -32    22   -12    421
   -40     4    -8     8    -8    -6   -36    46   -15    12    -65
    32    18   -32    -6   -11   -44     0     0     0     0   -337
   -32   -46   -24   -24    32    37    -8    42     0     3    124
   ];

[B1,P] = forwElimPPperm(A3);
P % permutation matrix
A = P*A3;
B = forwElim(A) % resulting echelon form


%% Question 9

clear all;
close all;
clc

% The code computes the approximations as provided by the Jacobi method and
% by the Gauss-Seidel method and presents the suppressed results in table form.

A = [...
    4 -1 0
    -1 8 -1
    0 -1 4
    ];
b = [48;12;24];
x0 = [1;1;1];
Nmax = 12;

% D is a diagonal matrix composed of the diagonal elements of A
D = diag(diag(A));

% L is composed of the elements below the main diagonal of A but with 
% opposite sign while U is composed of the elements above the main diagonal
% of A also with opposite sign
L = -1*tril(A,-1);
U = -1*triu(A,1);

% Jacobi method:

% Compute the T and c in Jacobi method
T1 = D\(L+U);
c1 = D\b;

% Compute the approximations and present them in a table
x_mat1 = genIter(T1,c1,x0,Nmax);
k = (0:Nmax)';
x1 = x_mat1(1,:)';
x2 = x_mat1(2,:)';
x3 = x_mat1(3,:)';
format long g
table(k,x1,x2,x3)


% Gauss-Seidel method:

% Compute the T and c in Gauss-Seidel method
T2 = (D-L)\U;
c2 = (D-L)\b;

% Compute the approximations and present them in a table
x_mat2 = genIter(T2,c2,x0,Nmax);
x1 = x_mat2(1,:)';
x2 = x_mat2(2,:)';
x3 = x_mat2(3,:)';
format long g
table(k,x1,x2,x3)


%% Question 10

% The code creats plots of l2 norm and l(inf) norm of the error v.s. k for 
% Jacobi method and Gauss-Seidel method each.

x = [13;4;7]; % exact solution x

% Jacobi method:

% Initialize two vectors norm_1 and norm_2 to store l2 and l(inf) norms
norm_1 = zeros(Nmax+1,1);
norm_2 = zeros(Nmax+1,1);

% Compute the error (denoted as error_1) between x and x(k) and its 
% l2 and l(inf) norm for each k

for i=1:Nmax+1 
    error_1 = x-x_mat1(:,i);
    norm_1(i) = norm(error_1,2);
    norm_2(i) = norm(error_1,Inf);
end

% Creat the figure
k = linspace(1,Nmax+1,Nmax+1);

subplot(2,1,1)
semilogy(k,norm_1,'b-')
hold on 
semilogy(k,norm_2,'r-')

title('The norm of error versus k for Jacobi method');
xlabel('k');
ylabel('norm of error');

legend('l2 norm','l(inf) norm');
hold off

% Fit polynomial to data
P1 = polyfit(k',norm_1,1);
P2 = polyfit(k',norm_2,1);

% Gauss-Seidel method:

% Initialize two vectors norm_3 and norm_4 to store l2 and l(inf) norms
norm_3 = zeros(Nmax+1,1);
norm_4 = zeros(Nmax+1,1);

% Compute the error (denoted as error_2) between x and x(k) and its 
% l2 and l(inf) norm for each k

for i=1:Nmax+1 
    error_1 = x-x_mat2(:,i);
    norm_3(i) = norm(error_1,2);
    norm_4(i) = norm(error_1,Inf);
end

% Creat the figure
subplot(2,1,2)
semilogy(k,norm_3,'b-')
hold on 
semilogy(k,norm_4,'r-')

title('The norm of error versus k for Gauss-Seidel method');
xlabel('k');
ylabel('norm of error');

legend('l2 norm','l(inf) norm');
hold off

% Fit polynomial to data
P3 = polyfit(k',norm_3,1);
P4 = polyfit(k',norm_4,1);

% From the plots, the slopes are -0.543, -0.473, -0.499 and -0.424 for
% l2 norm and l(inf) norm v.s. k for Jacobi method and Gauss-Seidel method 
% which can be obtained from P1, P2, P3 and P4.


%% Question 12

% The code creats the approximations provided by the Jacobi method and
% by the Gauss-Seidel method using genIterTol.m and presents the suppressed
% results in table form.

Nmax = 40;

tol = 1e-10;

% Jacobi method:
[x_mat3,Niter1,resNorm_vec1] = genIterTol(T1,c1,x0,Nmax,tol,A,b);

% Creat the table
format long g
k = (1:Niter1)';
table(k,resNorm_vec1)

% Gauss-Seidel method:
[x_mat4,Niter2,resNorm_vec2] = genIterTol(T2,c2,x0,Nmax,tol,A,b);

% Creat the table
format long g
k = (1:Niter2)';
table(k,resNorm_vec2)


%% Question 13

clear all;
close all;
clc

% The code in this section is separated into three parts which aim to 
% investigate a specific linear system Ax = b by Jacobi method. δ in the 
% question is replaced by p in the code.

% A and b are produced by function produce.m.

% Part 1: In this part, the code computes the T and c of the linear system
% (which is similiar to the one in Question 9). Here an additional function
% is built to produce the T and c for any linear system under the Jacobi
% method.

[A,b] = produce(0.5,7);
[T,c] = tAndc(A,b);


% Part 2: In this part, the code computes a figure of the total number of 
% iterations versus n which tends to infinity (100 for practice in this code).

clear all;
close all;
clc

tol = 1e-3;
Niter_vec = zeros(20,1);
Nmax = 50000;

for i = 1:20
    x0 = zeros(5*i,1);
    [A,b] = produce(0,5*i);
    [T,c] = tAndc(A,b);
    [x_mat,Niter,resNorm_vec] = genIterTol(T,c,x0,Nmax,tol,A,b);
    Niter_vec(i,1) = Niter;
end

% Creat the figure
n = linspace(1,100,20);
subplot(2,1,1)
plot(n,Niter_vec,'bo');
title('Number of iterations versus n');
xlabel('n');
ylabel('number of iterations');

% From the plot, we can see that the number of iterations increases to
% infinity as n goes to infinity.


% Part 3: In this part, the code computes figures of the total number of 
% iterations versus p for different Ns which tend to infinity.

% Creat the figures

Niter_vec = zeros(10,5);
Nmax = 50000;
p = linspace(1,10,10);
tol = 1e-3;
subplot(2,1,2)

for j = 1:5
    n = 100*j;
    for i = 1:10
        x0 = zeros(n,1);
        [A,b] = produce(i,n);
        [T,c] = tAndc(A,b);
        [x_mat,Niter,resNorm_vec] = genIterTol(T,c,x0,Nmax,tol,A,b);
        Niter_vec(i,j) = Niter;
    end
    
    plot(p,Niter_vec(:,j));
    hold on
end

title('Number of iterations versus p for several Ns');
xlabel('p');
ylabel('number of iterations');
legend('N1=100','N2=200','N3=300','N4=400','N5=500');

% From the plot, we can see that the number of iterations decreases to
% zero as p increases and this behaviour doesn't vary much as n tens to
% infinity.