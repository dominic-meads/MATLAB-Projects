%% Tutorial 1 - Sampling and Aliasing
% by Dominic Meads

clc;
clear all;
close all;
%% Example 1: High Sampling
% In this example we are going to plot a sinusoidal 
% signal and study the effect of sampling frequency

fmax = 3;               % max signal frequency
fs = 20*fmax;           % sampling frequency
Ts = 1/fs;              % sampling period
n = 0:1*fs-1;           % sample number
t = n*Ts;               % sampling instances
x = 10*cos(2*pi*3*t);

figure('Color', [1 1 1]);
plot(x);
xlabel('Sample Index');
ylabel('Signal');

figure('Color', [1 1 1]);
plot(t,x);
xlabel('Time (s)');
ylabel('Signal');

%% Example 2: Effects of Sampling
% Effects of reducing the sampling frequency

figure('Color', [1 1 1]);
h = plot(t,x);
set(h,'Linewidth',2);
set(h,'Color',[0.6 0.6 1]);
xlabel('Time (s)');
ylabel('Signal');
hold on;
h = plot(t,x,'.');
set(h,'Markersize',18);
set(h,'Color',[1 0 0]);

%% Example 3
% Effects of Aliasing

fmax = 1.5;               
fs = 20*fmax;           
Ts = 1/fs;              
n = 0:5*fs-1;           
t = n*Ts;               

x = 5 + 3*cos(2*pi*0.5*t) + 2*cos(1*pi*2*t) + 1*cos(2*pi*1.5*t);

figure('Color', [1 1 1]);
h = plot(t,x);
set(h,'Linewidth',2);
set(h,'Color',[0.6 0.6 1]);
xlabel('Time (s)');
ylabel('Signal');
hold on;

% undersample the signal, but LPF is applied first. 
fmax = 1.5;               
fs = 1.5;           
Ts = 1/fs;              
n = 0:5*fs-1;           
t = n*Ts;               
% all freuency components above 0.5*fs are removed
x1 = 5 + 3*cos(2*pi*0.5*t);

h = plot(t,x1,'.-');
set(h,'Markersize',18);
set(h,'Color',[1 0 0]);

% Sample signal at fs < 2*fmax
fmax = 1.5;               
fs = 1.5;           
Ts = 1/fs;              
n = 0:5*fs-1;           
t = n*Ts;               

x2 = 5 + 3*cos(2*pi*0.5*t) + 2*cos(1*pi*2*t) + 1*cos(2*pi*1.5*t);
figure('Color', [1 1 1]);
h = plot(t,x2);
set(h,'Linewidth',2);
set(h,'Color',[0 1 0]);
xlabel('Time (s)');
ylabel('Signal');


figure('Color', [1 1 1]);
h = plot(t,x2,'.-');
set(h,'Markersize',18);
set(h,'Color',[0 1 0]);
hold on;

% if we reconstruct the undersampled signal, it is not equal to x1. There
% is aliasing, which makes the reconstructed signal equal to x3.
fmax = 1.5;               
fs = 20*fmax;           
Ts = 1/fs;              
n = 0:5*fs-1;           
t = n*Ts;   
% reconstructed signal (some aliases add, some subtract, but overall a
% different signal than x1
x3 = 6 + 5*cos(2*pi*0.5*t);

h = plot(t,x3);
set(h,'Linewidth',2);
set(h,'Color',[0.6 0.6 1]);
xlabel('Time (s)');
ylabel('Signal');

%% Example 4: Real Signal
% Example with real world signal


% signal unavailable (no link)


% load icp.txt;
% figure('Color',[1 1 1]);
% icps = icp(10^4:10^4+1250); % segment
% fs = 125;
% t = (0:length(icps)-1)/fs;
% plot(t,icps);
% xlabel('Time (s)');
% ylabel('ICP (mmHg)');



