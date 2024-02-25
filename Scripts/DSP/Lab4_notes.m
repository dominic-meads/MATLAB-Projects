%% Lab 4 notes
% Dominic Meads

close all;
clear;
clc;

%% filter example
n = 1000; % order of filter
wn = 0.5; % normalized cuttoff frequency
b = fir1(n, wn, 'low');

stem(b); % filter coefficients (impulse response)
freqz(b); % frequency response


%% filter characteristics
load ('NoisyECG.mat');
fs = 500;
x = noisyECG_withTrend;

y = filter(b,1,x);


fs = 500;
f1 = 0;
f2 = 100;
n = 1000;

Wn1 = f1/(fs/2);
Wn2 = f2/(fs/2);

b = fir1(n,[0.4 0.5], 'high'); % bandpass


x1 = noisyECG_withTrend;
figure('Color', [1 1 1]);
t1 = (0:length(x1)-1)/fs;
noisyh = plot(t1,x1);
xlim([200 210]);
xlabel('Time (s)');
hold on;

y1 = filter(b,1,x1); % fixes time shift
h = plot(t1,y1);
xlim([200 210]);
set(h,'Color',[1 0.1 0.1]);
set(h,'Linewidth',2);

%%
stem(b);

title('Filter Coefficients');
xlabel('Sample Number');
ylabel('Coefficient Value');

[h,w] = freqz(b);
f = w/(2*pi)*fs;
plot(f,abs(h));
title('Frequency Response');
xlabel('Hz');
ylabel('Magnitude');

%% 
% smoothing use low pass,
% get rid of dc use high pass




