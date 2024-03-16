%% Tutorial - Event Detection
%
% Notes: if signal is more impulsive (dervative is high), use method in ECG
% QRS detection video. If signal is more rounded, use other things

close all
clear
clc
%% Objective
% Create a basic beat detection algorithm for cardiovascular pressure
% signals. 

%% Load signal

load ICP.mat

fs = 125;
icp = icp1;
t = (0:length(icp)-1)/fs;
figure('Color',[1 1 1]);
plot(t,icp);
xlabel('Time (s)');
ylabel('ICP (mmHG)');

% icp = icp2;
% t = (0:length(icp)-1)/fs;
% figure('Color',[1 1 1]);
% plot(t,icp);
% xlabel('Time (s)');
% ylabel('ICP (mmHG)');

%% Step 1 - Remove Low Frequency DC trend
% Apply a Highpass filter

fchp = 0.5; % respiration is a 0.5 hz or higher (record is 26 bpm low or 0.416 Hz)
Wp = fchp/(fs/2);
[B,A] = ellip(6,0.5,20,Wp,'high');

figure('Color',[1 1 1]);
freqz(B,A,2^10,fs);
figure('Color',[1 1 1]);
zplane(B,A);

icphp = filtfilt(B,A,icp);

figure('Color',[1 1 1]);
h = plot(t,icp);
hold on;
h = plot(t,icphp);

%% Step 2: LP Filter for Smoothing

fclp = 3;
Wp = fclp/(fs/2);
[B,A] = ellip(4,0.5,20,Wp);

figure('Color',[1 1 1]);
freqz(B,A,2^10,fs);
figure('Color',[1 1 1]);
zplane(B,A);

icplp = filtfilt(B,A,icphp);

figure('Color',[1 1 1]);
h = plot(t,icp);
hold on;
h = plot(t,icphp);
set(h,'Color',[0.5 1 0.5])
hold on;
h = plot(t,icplp);
set(h,'Color',[1 0.5 0.5])
hold on;

%% Step 3 Detection on Filtered Signal
x = icplp;
L = length(icplp);
id = find((x(1:L-2)<x(2:L-1))&(x(2:L-1)>x(3:L)))+1;
id2 = find((x(1:L-2)<=x(2:L-1))&(x(2:L-1)>=x(3:L)))+1;

h = plot(id/fs,x(id),'r.');
set(h,'MarkerSize',18);

figure('Color',[1 1 1]);
h = plot(t,icp,id/fs,icp(id),'r.');
set(h,'MarkerSize',18);

% you can overdetect, find all points, classify, then reject
figure('Color',[1 1 1]);
h = plot(t,icplp,id2/fs,icplp(id2),'r.');
set(h,'MarkerSize',18);
title('Overdetect');

%% Step 4 detection of first peak
% another option (1st option above)

% find max 
x = icplp;
L = length(icplp);
idmax = find((x(1:L-2)<x(2:L-1))&(x(2:L-1)>x(3:L)))+1;

% find min 
idmin = find((x(1:L-2)>x(2:L-1))&(x(2:L-1)<x(3:L)))+1;

figure('Color',[1 1 1]);
h = plot(t,icphp,t,icplp,idmax/fs,icphp(idmax),'r.',idmin/fs,icphp(idmin),'g.');
set(h,'MarkerSize',18);

% then take the interval between min and max of fundamental wave, apply
% interval to icp signal, and find the maximum in that interval (should
% give location of first peak

% Cardiac freqeuncy is fundamental freqeuncy





