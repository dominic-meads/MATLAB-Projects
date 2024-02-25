%% Tutorial 1 FIR Filtering
% by Dominic Meads
close all;
clear;
clc;
%% Concepts
% Moving average filters
% FIR Filters
% filter frequency response
%

%% Load Test Data
% Load ICP signal

load ICP.mat;
fs = 125;
x = icp1;

%% Plotting a real signal
% Plotting a signal with known fs

figure('Color',[1 1 1]);
t = (0:length(x)-1)/fs;
plot(t,x);
xlabel('Time(s)');
ylabel('ICP (mmHG)');
hold on;

%% Moving average filter 1
% smooth quantization noise
% 5-point window

M = 5;
B1 = 1/M*ones(M,1);

y1 = filtfilt(B1,1,x);
% filter is causal (depends on inputs x(present) and x(previous). It causes delay
% filtfilt is non-causal, and will not have time delay.
h = plot(t,y1);
set(h,'Color',[0.5 0.5 1]);
set(h,'Linewidth',2);
hold on;

%% Moving average filter 2
% smooth quantization noise to single sinusoid with large window

M = 30;
B2 = 1/M*ones(M,1);

y2 = filter(B2,1,x);
h = plot(t,y2);
set(h,'Color',[0.3 0.3 1]);
set(h,'Linewidth',2);

%% Frequency Response
% Comparing the frequnecy response
% A larger window will decrease cuttoff frequency

figure('Color',[1 1 1]);
freqz(B1,1,2^10,fs);
figure('Color',[1 1 1]);
freqz(B2,1,2^10,fs);

%% Design of an FIR Filter
% Using FIR1
fc = 8;
Wc = fc/(fs/2);
B3 = fir1(30,Wc);

y3 = filtfilt(B3,1,x);
h = plot(t,y3);
set(h,'Color',[1 0.5 0.5]);
set(h,'Linewidth',2);

%% Frequency Response of FIR filter
% Using freqz

figure('Color',[1 1 1]);
freqz(B3,1,2^10,fs);

%% Coefficients of filter
% AKA Impulse response
stem(B3);


%% IIR Filters
% Design and implementation of an elliptic filter

fc = 3;
Wc = fc/(fs/2);
[B,A] = ellip(8,0.5,20,Wc);
freqz(B,A,2^10,fs);

figure('Color',[1 1 1]);
y4 = filtfilt(B,A,x);
h = plot(t,y4);
set(h,'Color',[1 0.5 0.5]);
set(h,'Linewidth',2);















