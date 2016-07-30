%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script to discretely evaluate reactions using the Euler Algorithm

%clc;
clear all;
%close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% declare geometry and dynamics
V = 2.5*10^4;
A = 4.4*10^3;
NA = 2.4*10^5;
NP = 9.8*10^4;
DA = 0.28;
DP = 0.15;
% declare reaction rates
koffA = 3.24*10^(-3); % /s
koffP = 7.19*10^(-3); % /s
konA = 6.29*10^(-3); % um/s
konP = 7.682*10^(-2); %um/s
kAP = 0.004; % um^2/s % measure is obtained from other groups
kPA = 0.006; % um^2/s % measure is obtained from other groups

a=1;
b=2;

t = 1;
tmax = 500;
dt = 0.05;
tsteps = tmax/dt;

% declare initial conditions
Am0 = (0.45*10^5)/A;      %0.05*NA;
Ac0 = NA-Am0*A;
Pm0 = (0.05*10^5)/A;      %0.45*NP;
Pc0 = NP-Pm0*A;

Am = zeros(1,tsteps);
Am(1) = Am0;

Ac = zeros(1,tsteps);
Ac(1) = Ac0;

Pm = zeros(1,tsteps);
Pm(1) = Pm0;

Pc = zeros(1,tsteps);
Pc(1) = Pc0;

% solve reaction equations
for t = 1:tsteps
    %Am(t+1) = Am(t)*(dt*(koffA+A/V+Pm(t).^a*kAP)+1)-NA/V*konA*dt;
    Am(t+1) = Am(t)-Am(t)*koffA*dt+konA*dt*((NA-Am(t)*A)/V)-Am(t)*Pm(t).^a*kAP*dt;
    Ac(t+1) = NA/A - Am(t+1);
    %Pm(t+1) = Pm(t)*(1-dt*(koffP+A/V*konP+Am(t).^b*kPA))+NP/V*konP*dt;
    Pm(t+1) = Pm(t)-Pm(t)*koffP*dt+konP*dt*((NP-Pm(t)*A)/V)-Pm(t)*Am(t).^b*kPA*dt;
    Pc(t+1) = NP/A - Pm(t+1);
end

tvector = [1:(tsteps+1)];
realtime = tvector*dt;

figure(1);
plot(realtime,(Am/A));

%figure(3);
%plot(realtime,Pm);

