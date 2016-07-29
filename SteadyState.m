% Steady state analysis

clear all;

V = 2.5*10^4;
O = 4.4*10^3;
NA = 2.4*10^5;
NP = 9.8*10^4;
DA = 0.28;
DP = 0.15;
koffA = 3.24*10^(-3);
koffP = 7.19*10^(-3);
konA = 6.29*10^(-3);
konP = 7.682*10^(-2);
kAP = 3.24*10^(-1);
kPA = 7.19*10^(-1);

%Am1 = ((NA.*konA)./(koffA+O/V+Pm.^a.*kAP));
%Am2 = ((((NP-Pm.*O)./V).*konP-Pm.*koffP)./kPA.*Pm).^(1./b);

a = 2;
Pm = linspace(0, 30, 100);
for b= 1            % b = (1:0.5:3)
    Am1 = ((NA./V.*konA)./(koffA+O/V+Pm.^(a).*kAP));
    Am2 = ((((NP-Pm.*O)./V).*konP-Pm.*koffP)./(kPA.*Pm)).^(1./b);
    plot(Pm, Am1)
    hold on
    plot(Pm, Am2)
    hold on    
end
