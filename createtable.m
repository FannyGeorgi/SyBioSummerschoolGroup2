data=importdata('ExperimentalData\T_PARs.txt');
T(:,1)=data';

data=importdata('ExperimentalData\X_PARs.txt');
S=data';

data=importdata('ExperimentalData\PAR2.txt');
P2=data;

data=importdata('ExperimentalData\PAR6.txt');
P6=data;
%Labelt='time[s]';
HeatMap(P6,'RowLabels',T)

