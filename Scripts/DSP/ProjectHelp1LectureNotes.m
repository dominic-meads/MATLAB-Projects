%% Event Detection 

% Objective is to implement a basic QRS detection algorithm
% The QRS is inspired on the Pan-Tompkins algorithm for QRS Detection

% Algorithm 
% Step 1 - Lowpass Filter (60 Hz, smoothing)
% Step 2 - Highpass Filter
% Step 3 - Differentiator (enhance QRS complex)
% Step 4 - Squaring Operator
% Step 5 - Moving Average Filter (integration)
% Step 6 - Decision Logic

close all
clear
clc
%% Load Noisy ECG
% load Noisy ECG signal sampled at 500 Hz
load ECGNoisy60Hz.mat

% Decimate, Pan-Tompkins algorithm designed for fs = 200 Hz
ecg = decimate(ecgn,2);
ecg = ecg-mean(ecg);
fs = fs/2;

%% Plot ECG signal

t = (0:length(ecg)-1)/fs;
figure('Color', [1 1 1]);
h = plot(t,ecg);
xlabel('Time (s)');
ylabel('ECG');
box off;
hold on;

%% Step 1
% Lowpass filter to Eliminate High frequency noise and 60 Hz noise

fclp = 22; % 11 in the book (work both for adults and children)
Wp = fclp/(fs/2);
[B,A] = ellip(8,0.5,20,Wp);
% freqz(B,A,2^10,fs);
% zplane(B,A);

ecglp = filtfilt(B,A,ecg);
h = plot(t,ecglp);
set(h,'Color',[0.5 0.5 1]);
set(h,'LineWidth',2);

%% Step 2:
% Highpass filter to eliminate low frequency noise (Baseline drift)

fchp = 5;
Wp = fchp/(fs/2);
[B,A] = ellip(8,0.5,20,Wp,'high');
% freqz(B,A,2^10,fs);
% zplane(B,A);

ecghp = filtfilt(B,A,ecglp);
h = plot(t,ecghp);
set(h,'Color',[1 0.5 0.5]);
set(h,'LineWidth',2);
title('Step 1 and 2: LP and HP');

%% Step 3: Differentiator
% derviative operators to "ehnance" QRS Feature

ecgd = diff([0; ecghp]); % append zero to beginning to make signals same length
figure('Color',[1 1 1]);
h = plot(t,ecgd);
set(h,'Color',[0.5 0.5 1]);
set(h,'LineWidth',2);
title('Step 3: Differentiator');

%% Step 4: Square Operator

ecgs = ecgd.^2; % append zero to beginning to make signals same length
figure('Color',[1 1 1]);
h = plot(t,ecgs);
set(h,'Color',[0.5 0.5 1]);
set(h,'LineWidth',2);
title('Step 4: Square Operator')

%% Step 5: Moving Average Filter
% Integrator

B = ones(30,1)/30;
ecgm = filtfilt(B,1,ecgs);
figure('Color',[1 1 1]);
h = plot(t,ecgm);
set(h,'Color',[0.5 0.5 1]);
set(h,'LineWidth',2);
title('Step 5: MA')


%% Step 6: Decision Logic
% Peak detection

x = ecgm;
nx = length(x);
id = find((x(1:nx-2)<x(2:nx-1))&(x(2:nx-1)>x(3:nx)))+1;
id2 = find(ecghp(id) > 1);
id2 = id(id2);
figure('Color',[1 1 1]);
h = plot(t,x,id/fs,x(id),'r.');
set(h,'MarkerSize',18);
title('Step 6: Decision Logic');

%% Step 7: Detection on Original ECG
%
% If point is not at top, peak detection is off by 1 sample or 1/25 = 4 ms
% Want dot to touch QRS complex

figure('Color',[1 1 1]);
h = plot(t,ecg,id/fs,ecg(id),'r.');
set(h,'MarkerSize',18);
set(h,'LineWidth', 2);
title('Step 7: Detection on Original Signal');
xlabel('Time (s)');
ylabel('ECG');
hold on;

h = plot(t,ecg,id2/fs,ecg(id2),'g.');

%% Step 8: Interbeat Interval Analysisfind overdetections

figure('Color',[1 1 1]);
ibi = diff(id);
h = plot(ibi);
set(h,'LineWidth', 2);
title('Step 8: Interbeats');
xlabel('sample');
ylabel('IBI');

% in graph, downward peaks mean small "interbeats" which is over detection.
% Heartrate doestn change that quickly. Improved by high fc in LP filter
% (not modifying signal as much)

figure('Color',[1 1 1]);
ibi2 = diff(id2);
h = plot(ibi2);
set(h,'LineWidth', 2);
title('Step 8: Interbeats');
xlabel('sample');
ylabel('IBI');







