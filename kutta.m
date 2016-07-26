#runge kutta
#dc/dt = -lambda*c     c(t = 0) =c0 

#assume lambda = 1

clc
close all
clear all
global lambda =1;
x0 = 1;
t0 = 0;
h = 0.3; #h is my dt
i = 1;
n = 100;
x(i) = x0;
while i<=n
  k1 = h*f1(t0,x0);
  k2 = h*f1(t0+h/2,x0+k1/2);
  k3 = h*f1(t0+h/2,x0+k2/2);
  k4 = h*f1(t0+h,x0+k3);
  km = (k1+2*k2+2*k3+k4)/6;
  t1 = t0 + h;
  x1 = x0 + km;
  i = i + 1;
  x(i) = x1;
  t(i) = t1;
  x0 = x1;
  t0 = t1;
end

plot(t,x)
grid on