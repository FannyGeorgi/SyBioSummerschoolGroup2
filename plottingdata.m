%%%% script to generate plots for experimental data of SysBioSummerSchool Dresden 2016

% import txt files containing timepoints, distances and intensities at distances over time
data=importdata('ExperimentalData/T_PARs.txt'); % mind that using Windows, \ is used instead of /
%T(:,1)=data';                                  % alternative way for creating array
T=data';
data=importdata('ExperimentalData/X_PARs.txt');
S=data';
data=importdata('ExperimentalData/PAR2.txt');
P2=data;
data=importdata('ExperimentalData/PAR6.txt');
P6=data;

% attempt to use HeatMap function
%Labelt='time[s]';
%HeatMap(P6,'RowLabels',T)

% attempt to use imagesc
colorscheme = othercolor('Blues9'); 
imagesc(P2)
colorbar;

colormap('hot');                                  % set colormap
imagesc(P6);                                      % draw image and scale colormap to values range
colorbar;                                         % show color scale

