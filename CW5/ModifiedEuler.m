function [t,w]=ModifiedEuler(a,b,N,y0,f)
%%This funtion computes approximation by Modified Euler's Method to solve 
%%IVP on the interval [a,b] with N+1 mesh points. y0 is the initial 
%%condition and f is the right handside of the ODE in the IVP. The output 
%%t is a vector with meshpoints and w a vector with approximations at these
%%meshpoints.

h=(b-a)/N;
w=zeros(N+1,1); % initialise
t=zeros(N+1,1);
w(1) = y0;
t(1)=a;
for i=1:N
  t(i+1) = t(i) + h;
  w(i+1) = w(i) + (h/2) * (f(t(i),w(i)) + f(t(i+1),w(i)+h*f(t(i),w(i))));
end