%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script to treat reactions (and diffusion) probabilistic using the Guiespie Algorithm

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
tmax = 500000;

% declare initial conditions
Am0 = (0.45*10^5);      %0.05*NA;
Ac0 = NA-Am0;
Pm0 = (0.05*10^5);      %0.45*NP;
Pc0 = NP-Pm0;

N0 = [Am0, Ac0, Pm0, Pc0];
Nt = N0';
Nt_arr = zeros(4,tmax+1);
Nt_arr(:,1) = Nt;
%Nt = 100*ones(1,4);

% declare a and b
a = 1;
b = 2;

% declare master vectors contains production (Vplus) and consumption (Vminus) during reactions (colums)
Vplus = [0, 1, 0, 0, 0, b; 1, 0, 0, 0, 1, 0; 0, 0, 0, 1, a, 0; 0, 0, 1, 0, 0, 1];
Vminus = [1, 0, 0, 0, 1, b; 0, 1, 0, 0, 0, 0; 0, 0, 1, 0, a, 1; 0, 0, 0, 1, 0, 0];

% declare reaction rates
koffA = 3.24*10^(-3); % /s
koffP = 7.19*10^(-3); % /s
konA = 6.29*10^(-3); % um/s
konP = 7.682*10^(-2); %um/s
kAP = 0.004; % um^2/s % measure is obtained from other groups
kPA = 0.006; % um^2/s % measure is obtained from other groups

k = [koffA, koffP, konA, konP, kAP, kPA];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate specific probability rates of reactions
c= zeros(1,6);
for mu=1:6
    c(mu) = k(mu).*factorial(Vminus(1, mu)).*k(mu).*factorial(Vminus(2, mu)).*k(mu).*factorial(Vminus(3, mu)).*k(mu).*factorial(Vminus(4, mu));
end
c = [c(1), c(2), c(3).*V./A, c(4).*V./A, c(5).*A.^(-a), c(6).*A.^(-b)];

% iterate over time steps
h = zeros(1,6);
ttime = zeros(1,tmax);
ttime_sum = zeros(1,tmax);
reaction = zeros(1,tmax);
for t = 1:tmax
     % calculate number of possible collisions
    for mu = 1:6
        h(mu) = nchoosek(Nt(1), Vminus(1, mu))*nchoosek(Nt(2), Vminus(2, mu))*nchoosek(Nt(3), Vminus(3, mu))*nchoosek(Nt(4), Vminus(4, mu));
    end
    
    % calculate propensity
    aprop = c.*h;
    adot = sum(aprop);
    
    % evaluate the probability for each reaction
    p = aprop./adot;
    
    % calculate randomized time step duration
    tau = -log(rand(1))/adot;
    
    % chose reaction that will happen in this time step
    check = p(1);
    reactionhappening = 2;
    dice = rand(1);
    while dice > check
        check = check + p(reactionhappening);
        %disp(check);
        reactionhappening = reactionhappening +1;
    end
    
    % let the reaction happen
    Nt = Nt + Vplus(:,reactionhappening-1) - Vminus(:,reactionhappening-1);
    Nt_arr(:,(t+1)) = Nt;
    
    % keep track of reacktions happening
    reaction(t) = reactionhappening-1;
    
    % save duration of time step
    ttime(t) = tau;
    if t == 1
        ttime_sum(t) =tau;
    else 
        ttime_sum(t) = ttime_sum(t-1) + tau;
    end
    
end

figure(1);
plot([0,ttime_sum],Nt_arr(1,:));

%figure(2);
%plot([0,ttime_sum],Nt_arr(2,:));

%figure(3);
%plot([0,ttime_sum],Nt_arr(3,:));

%figure(4);
%plot([0,ttime_sum],Nt_arr(4,:));    
