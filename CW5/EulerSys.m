function [t,w1,w2]=EulerSys(T,N,theta,alpha,delta)
%%This funtion computes Euler's approximation to the solution of a system of
%%two first-order IVPs on the interval [0,T] with N+1 mesh points. alpha 
%%and delta are two initial conditions and theta is a parameter. The output 
%%t is a vector with meshpoints. w1 and w2 are two vectors with approximations
%%to the system satisfied by (u1,u2) at these meshpoints.

h=T/N;
w1=zeros(N+1,1); % initialise
w2=zeros(N+1,1);
t=zeros(N+1,1);
w1(1) = alpha;
w2(1) = delta;
t(1)=0;
% The system of ODE's
f1 = @(u2) u2;
f2 = @(u1) -theta^2*u1;
for i=1:N
  w1(i+1) = w1(i) + h * f1(w2(i));
  w2(i+1) = w2(i) + h * f2(w1(i));
  t(i+1) = t(i) + h;
end



