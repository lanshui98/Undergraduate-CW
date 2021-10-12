function [t,w]=Euler(a,b,N,y0,f)
%%This funtion computes Euler approximation to solve IVP on the interval
%%[a,b] with N+1 mesh points. y0 is the initial condition and f is the
%%right handside of the ODE in the IVP. The output t is a vector with 
%%meshpoints and w a vector with Euler approximations.

h=(b-a)/N;
w=zeros(N+1,1);
t=zeros(N+1,1);
w(1) = y0;
t(1)=a;
for i=1:N
  w(i+1) = w(i) + h * f(t(i),w(i));
  t(i+1) = t(i) + h;
end