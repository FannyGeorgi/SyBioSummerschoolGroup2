clc; clear all; close all;

min=0;
max=100;
dt=0.1;
c0=50;
l=1;

N=(max-min)/dt;
t1=linspace(min,max,N);
t2=linspace(min,max,1000);

c=zeros(1,N);
c(1)=c0;

for i=1:N-1
    c(i+1)=c(i)*(1-l*dt);
end

th=c0*exp(-t2);
plot(t1,c)
hold on
plot(t2,th,'r')
axis([0 8 0 100]);