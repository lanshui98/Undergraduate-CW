function [t,w] = RuKuMeth(a,b,N,y0,f)
%%This funtion computes approximation by Runge–Kutta Order Four method to 
%%solve IVP on the interval [a,b] with N+1 mesh points. y0 is the initial 
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
  K1 = h * f(t(i),w(i));
  K2 = h * f(t(i)+h/2,w(i)+K1/2);
  K3 = h * f(t(i)+h/2,w(i)+K2/2);
  K4 = h * f(t(i)+h,w(i)+K3);
  w(i+1) = w(i) + (K1 + 2*K2 + 2*K3 + K4)/6;
end






