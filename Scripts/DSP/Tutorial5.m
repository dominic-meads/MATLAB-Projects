%% Tutorial - Spectral analysis on real signals
% by Dominic Meads
%
% When analyzing real signals, decimate, remove DC trend, and window
close all
clear
clc
%% Load Signal
load ICP.mat
x = icp1;

%% Plot signal
figure('Color',[1 1 1]);
fs = 125;
t = (0:length(x)-1)/fs;
h = plot(t,x);
xlabel('Time (s)');
ylabel('ICP (mmHG)');

%% Spectrum

N = 2^12;
x = x.*blackman(length(x));
X = abs(fft(x,N));
X = X(1:end/2);
f = linspace(0,fs/2,N/2);

figure('Color',[1 1 1]);
h = plot(f,X);
xlabel('Frequency (Hz)');
ylabel('Spectrum');

% before running FFT, decimate and remove DC component (high pass filter)
x = decimate(x,5); % decimate by factor of 5
fs = fix(fs/5);

% quick tricks to remove DC:
% x = x-mean(x); % remove average (DC value of signal) ONLY FOR SMALL SECTIONS OF SIGNAL
% x = diff(x); % derivative, removes low freqeuency but amplifies high frequnecy

% high pass filter
% easier to make a stable filter when max frequency is 10*fcuttoff, we
% decimated to make fmax = 12.5
Wp = 0.3/(fs/2);
[B,A] = ellip(6,0.5,20,Wp,'high');
% figure('Color',[1 1 1]);
% freqz(B,A,2^10,fs);
x = filtfilt(B,A,x);

% z plane
% figure('Color',[1 1 1]);
% zplane(B,A);

X = abs(fft(x.*blackman(length(x)),N));
X = X(1:end/2);
f = linspace(0,fs/2,N/2);

figure('Color',[1 1 1]);
h = plot(f,X);
xlabel('Frequency (Hz)');
ylabel('Spectrum');
title('Decimated and DC component removed');
