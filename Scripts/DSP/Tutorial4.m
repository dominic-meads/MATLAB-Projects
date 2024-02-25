%% Tutorial on Spectrum analysis
% by Dominic Meads
close all
clear 
clc

%% Create signal
% x(t) = cos(2*pi*5*t)

fmax = 5;
fs = 5*fmax;
Ts = 1/fs;
n = 0:1*fs;
t = n*Ts;
x = cos(2*pi*5*t);

figure('Color',[1 1 1]);
h =  plot(t,x);
xlabel('Time (s)');
ylabel('x(t)');

%% Practice using fft
% one-sided spectrum

% Magnitude (usually dont need more than 2^12 points, 2^10 most common)
N = 2^12; % must use power of 2 for fft (otherwise uses DFT)
X = abs(fft(x,N));

% make one sided
X = X(1:end/2);

% create freqeuncy vector to plot
% create N/2 pooints between min freq and max freq (nyquist)
f = linspace(0,fs/2,N/2);
% addtional way to create frequency vector
% f = fs/N*(0:N/2-1);

figure('Color',[1 1 1]);
plot(f,X);
% plot(f,20*log(X));  % to see signals with very different ratios
xticks(0:5:fs/2);
txt = {'\leftarrow Note side lobes due to', 'rectangular window'};
text(8.5,2.5,txt,'FontSize',10);

%% Computational Resolution
% Number of points in the FFT

N = 2^12; %% Keep number of points higher for better COMPUTATIONAL resolution
X = abs(fft(x,N));

X = X(1:end/2);

f = linspace(0,fs/2,N/2);

figure('Color',[1 1 1]);
subplot(3,1,1);
plot(f,X);
xticks(0:5:fs/2);
title('High Computational Resolution (N = 2^{12})');
hold on; 

% Note, to "Zoom in" while using the same number of points,
% you can decrease the sampling frequency (decrease max
% x axis in fft plot)

%% Frequency Resolution
% How long we are looking at the signal

n = 0:0.3*fs; % Rectangular window is smaller
t = n*Ts;
x = cos(2*pi*5*t);
N = 2^12; 
X = abs(fft(x,N));
X = X(1:end/2);
f = linspace(0,fs/2,N/2);
subplot(3,1,2);
plot(f,X);
xticks(0:5:fs/2);
title('Low Frequency Resolution (Time-domain signal only plotted for 0.3 sec)');

% look for longer

n = 0:10*fs; % Rectangular window is much larger (10 sec)
t = n*Ts;
x = cos(2*pi*5*t);
N = 2^12; 
X = abs(fft(x,N));
X = X(1:end/2);
f = linspace(0,fs/2,N/2);
subplot(3,1,3);
plot(f,X);
xticks(0:5:fs/2);
title('High Frequency Resolution (Time-domain signal plotted for 10 sec)');

%% Windowing 
% Different windows
% * Blackman (lose a little frequency resolution, but eliminate side lobes)
% * Hamming
% * Hanning
% * rectangular
% * etc...

w = blackman(length(x));
figure('Color',[1 1 1]);
h = plot(w);
title('Blackman Window');

% window the signal
xw = x.*w';
figure('Color',[1 1 1]);
h = plot(xw);
title('Signal Multiplied By Blackman Window');

% comput fft
N = 2^12; 
Xw = abs(fft(xw,N));
Xw = Xw(1:end/2);
f = linspace(0,fs/2,N/2);
figure('Color',[1 1 1]);
plot(f,Xw);
xticks(0:5:fs/2);
txt = {'\leftarrow No side lobes', 'Blackman window'};
text(6.5,10,txt,'FontSize',14);

%% FFT of different window

win = ["Blackman" "Hanning" "Hamming" "Rectangle" "Triangle"];
w = zeros(5,length(x));
w(1,:) = blackman(length(x));
w(2,:) = hann(length(x));
w(3,:) = hamming(length(x));
w(4,:) = rectwin(length(x));
w(5,:) = triang(length(x));

figure('Color',[1 1 1]);
for k = 1:5
    subplot(3,2,k);
    plot(w(k,:));
    title(strcat(win(k),' Window'));
end


% ffts of windows
N = 2^12;
figure('Color',[1 1 1]);
X = zeros(5,N);

figure('Color',[1 1 1]);
for k = 1:5
    subplot(3,2,k);
    X(k,:) = abs(fft(w(k,:),N));
    plot(X(k,:));
    title(strcat('{FFT of }',win(k),' Window'));
end

%% Frequency Response

% Moving average filter Impulse response
B = ones(10,1)/10;

% compute fft and compare to frequnecy response (they are the same)
N = 2^12; 
X = abs(fft(B,N));
X = X(1:end/2);
f = linspace(0,fs/2,N/2);
figure('Color',[1 1 1]);
plot(f,20*log(X));
figure('Color',[1 1 1]);
freqz(B,1,2^10,fs);









