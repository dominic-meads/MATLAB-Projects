%% FFT magnitude notes

close all
clear 
clc

%% create sinusoid 
% DC component = 3
% magnitude 10 @ f = 10 Hz

fs = 100;
Ts = 1/fs;
n = 0:1*fs-1;
t = n*Ts;
x = 3 + 10*cos(2*pi*10*t);

figure('Color',[1 1 1]);
h = plot(t,x);
title('x(t)');
ylabel('Amplitude');
xlabel('Time (s)');

%% plot single sided spectrum 
% Expected magnitudes:
% * 10 @ f = 10 Hz
% * 3 @ f = 0 Hz

N = 2^10;
X1 = fft(x,N);  
X1 = X1/length(x);              % normalize by dividing by length of input signal
X1 = X1(1:end/2);               % single-sided
X1(2:end) = 2*X1(2:end);        % fix halved magnitude (exclude DC)
f1 = linspace(0,fs/2,N/2);
figure('Color',[1 1 1]);
h = plot(f1,abs(X1));           % Plot magnitude
title('X(jw) Single-Sided Spectrum');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

%% plot double sided spectrum 
% Expected magnitudes:
% * 5 @ f = 10 Hz
% * 5 @ f = -10 Hz
% * 3 @ f = 0 Hz

N = 2^10;
X2 = fft(x,N);  
X2 = X2/length(x);                % normalize by dividing by length of input signal
X2 = [X2(N/2+2:end) X2(1:N/2+1)]; % split data in the middle, then flip each section
f2 = linspace(-fs/2,fs/2,N);      % update linspace to include negative frequencies
figure('Color',[1 1 1]);
h = plot(f2,abs(X2));             % Plot magnitude
title('X(jw) Double-Sided Spectrum');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

%% plot single sided spectrum With Decimation 
% Same method. 
% Expected magnitudes:
% * 10 @ f = 10 Hz
% * 3 @ f = 0 Hz

x = decimate(x,2);
fs = fix(fs/2);

N = 2^10;
X1 = fft(x,N);  
X1 = X1/length(x);              % normalize by dividing by length of input signal
X1 = X1(1:end/2);               % single-sided
X1(2:end) = 2*X1(2:end);        % fix halved magnitude (exclude DC)
f1 = linspace(0,fs/2,N/2);
figure('Color',[1 1 1]);
h = plot(f1,abs(X1));           % Plot magnitude
title('X(jw) Single-Sided Spectrum');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
