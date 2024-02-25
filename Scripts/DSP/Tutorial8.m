%% Module 5
% Non stationary analysis
close all
clear 
clc
%% Simple Example of the Spectrogram
% Spectrogram of a sinusoidal signal 

% Sinusoidal signal
close all;
fs = 125;
Ts = 1/fs;
tn = 0:Ts:10;
f = 2;
x = 2*cos(2*pi*f*tn)+1*cos(2*pi*2*f*tn+pi/2);
figure('Color',[1 1 1]);
h = plot(tn,x);
xlabel('Time (s)');
ylabel('Amplitude');

% Spectrogram 
[y,f,t,p] = spectrogram(x,125,[],2^12,fs,'yaxis');
figure('Color',[1 1 1]);
imagesc(t,f,p);
set(gca, 'YDir', 'normal');
xlabel('Time (s)');
ylabel('Frequency (Hz)');
colorbar

%% Spectrogram of decimated signal
xd = decimate(x,10);
fsd = fix(fs/10);
[y,f,t,p] = spectrogram(xd,12.5,[],2^12,fsd,'yaxis');
figure('Color',[1 1 1]);
imagesc(t,f,p);
set(gca, 'YDir', 'normal');
xlabel('Time (s)');
ylabel('Frequency (Hz)');
colorbar

figure('Color',[1 1 1]);
h = plot(f,mean(p'));
xlabel('Frequency (Hz)');
ylabel('PSD');

%% Spectrogram of decimated signal (increased window size)
xd = decimate(x,10);
fsd = fix(fs/10);
[y,f,t,p] = spectrogram(xd,25,[],2^12,fsd,'yaxis');
figure('Color',[1 1 1]);
imagesc(t,f,p);
set(gca, 'YDir', 'normal');
xlabel('Time (s)');
ylabel('Frequency (Hz)');
colorbar

figure('Color',[1 1 1]);
h = plot(f,mean(p'));
xlabel('Frequency (Hz)');
ylabel('PSD');

%% Spectrogram of decimated signal 3D plot
xd = decimate(x,10);
fsd = fix(fs/10);
[y,f,t,p] = spectrogram(xd,25,[],2^12,fsd,'yaxis');
figure('Color',[1 1 1]);
meshc(t,f,p);
view(0,90);
set(gca, 'YDir', 'normal');
xlabel('Time (s)');
ylabel('Frequency (Hz)');
colorbar

%% Spectrogram of decimated signal 3D in dB
xd = decimate(x,10);
fsd = fix(fs/10);
[y,f,t,p] = spectrogram(xd,25,[],2^12,fsd,'yaxis');
figure('Color',[1 1 1]);
meshc(t,f,10*log(p));
view(0,90);
set(gca, 'YDir', 'normal');
xlabel('Time (s)');
ylabel('Frequency (Hz)');
colorbar

%% Spectrogram of real biomedical signal
% ICP signal

load ICP.mat
fs = 125;
x = icp1;
t = (0:length(x)-1)/fs;
figure('Color',[1 1 1]);
h = plot(t,x);
title('ICP1 signal');
xlabel('Time(s)');
ylabel('ICP (mmHg)');

% step 1 decimate signal
dec = 10;
xd = decimate(x,dec);
fsd = fix(fs/dec);

% step 2 remove local trend
Wp = 0.1/(fsd/2);
[B,A] = ellip(6,0.5,20,Wp,'high');
xd = filtfilt(B,A,xd);

% Spectrogram
[y,f,t,p] = spectrogram(xd(1:2000),25,[],2^12,fsd,'yaxis');
figure('Color',[1 1 1]);
% meshc(t,f,10*log(p));
meshc(t,f,p);
view(0,90);
set(gca, 'YDir', 'normal');
title('Spectrogram of ICP1');
xlabel('Time (s)');
ylabel('Frequency (Hz)');
colorbar


% Movie
figure('Color',[1 1 1]);
for k = 1:length(t)
    plot(f,p(:,k));
    pause(0.1);
end


% average (spectrogram function of whole signal
figure('Color',[1 1 1]);
plot(f,mean(p'));




