%% Lab 5 
% by Dominic Meads

close all
clear
clc

%% TODO
%
% ALL SECTIONS
% FIX FFT MAGNITUDE (DIVIDE BY LENGTH OF SIGNAL)
%


%% 1). FFT of Signal With Different Window Lengths
%
%  x(t) = 3sin(2pi*0.5t) + 2sin(2pi*2.5t) 
%
% Increased window length results in greater frequency resolution. The
% peaks in the spectrum are sharper and more closley resemble impulses,
% however, we do see some sidelobes. When the window is small, we do not see
% many sidelobes, but the peaks are wide (if two signals were close
% together in frequency, it would be difficult to see).

% ------------------ Plot signal --------------------------------------
fmax = 2.5;
fs = 5*fmax;
Ts = 1/fs;
n = 0:2*fs;
t = n*Ts;
x = 3*sin(2*pi*0.5*t) + 2*sin(2*pi*2.5*t);
figure('Color',[1 1 1]);
sgtitle('Spectrum Using Rectanglular Window');
subplot(2,2,1);
h = plot(t,x);
title('x(t)=3sin(2\pi0.5t)+2sin(2\pi2.5t)');
xlabel('Time (s)');
ylabel('Amplitude');

% ------------------- Calculate FFT at different window sizes --------
win_size = [1 20 100];
N = 2^12;
f = linspace(0,fs/2,N/2);

for k = 1:numel(win_size)
    n = 0:win_size(k)*fs;
    t = n*Ts;
    x = 3*sin(2*pi*0.5*t) + 2*sin(2*pi*2.5*t);
    X = abs(fft(x,N));
    X = X(1:end/2);   % plot single-sided spectrum
    subplot(2,2,1+k);
    h = plot(f,X);
    title(strcat('{Spectrum x(t) (win length = }', ...
        num2str(win_size(k)),')'));
    xlabel('Frequency (Hz)');
    ylabel('Spectrum');
end

%% 2). Comparing Windowed FFTs of a Signal With Different Window Lengths
%
% The triangular window has the sharpest peaks when window length is high,
% however, it has the 2nd most sidelobes (rectangular has the most). You
% can still make out the freqeuncy components in the 1-second length
% Blackman window, but it also has one small sidelobe. Additionally, The peaks are
% not as sharp as in the triangular example. Finally, the Blackmanharris
% window has no side lobes, but in short window lengths you cannot discern
% induvidual frequency components that are close together. 

win_type = {'Triangular','Blackman','Blackmanharris'};
win_fun = {'triang','blackman','blackmanharris'};
w = cellfun(@str2func,win_fun,'uni',0);

for i = 1:numel(w)      % Calculates ffts for window types
    figure('Color',[1 1 1]);
    sgtitle(strcat('Spectrum of x(t) (', win_type(i), ' window)'));
    for k = 1:numel(win_size)      % calculate ffts for window lengths
        n = 0:win_size(k)*fs;
        t = n*Ts;
        x = 3*sin(2*pi*0.5*t) + 2*sin(2*pi*2.5*t);
        x = x'.*(w{i}(numel(x)));    % Call function
        X = abs(fft(x,N));
        X = X(1:end/2);             % plot single-sided spectrum
        subplot(3,1,k);
        h = plot(f,X);
        title(strcat('{Spectrum of x(t) (win length = }', ...
            num2str(win_size(k)),')'));
        xlabel('Frequency (Hz)');
        ylabel('Spectrum');
    end
end

%% 3). Generating A Sinusoid and Visualizing Window Effects
%
% It is difficult to see the sidelobes (effect of rectangular window) in the figure. 
% They could be more apparent if the signal was decimated. This would 
% lower the maximum frequency on the x-axis in the FFT plot, effectivley
% "zooming in" without losing resolution. 

% ------------------ Plot signal --------------------------------------
fmax = 8;
fs = 64*fmax;
Ts = 1/fs;
n = 0:1*fs-1;  % 512 samples (rectangular window)
t = n*Ts;
x = sin(2*pi*8*t);
figure('Color',[1 1 1]);
h = plot(t,x);
title('x(t)=sin(2\pi8t)');
xlabel('Time (s)');
ylabel('Amplitude');

% -------------------- FFT ------------------------------------------

N = 2^12;
f = linspace(0,fs/2,N/2);
X = abs(fft(x,N));
X = X(1:end/2);
figure('Color',[1 1 1]);
h = plot(f,X);
title('FFT of x(t)=sin(2\pi8t)');
xlabel('Frequency (Hz)');
ylabel('Spectrum');

%% 4). Problem 1 with different N-point FFTs

% ------------------ Plot signal --------------------------------------
fmax = 2.5;
fs = 5*fmax;
Ts = 1/fs;
n = 0:2*fs;
t = n*Ts;
x = 3*sin(2*pi*0.5*t) + 2*sin(2*pi*2.5*t);

% ------------------- Calculate different N-point FFTs ------------
win_size = [1 20 100];
N = [2^8 2^10 2^12];

for i = 1:numel(N)      % Calculates ffts for different N
    figure('Color',[1 1 1]);
    sgtitle(strcat('FFT using N=', num2str(N(i)), ' points)'));
    f = linspace(0,fs/2,N(i)/2);
    for k = 1:numel(win_size)      % calculate ffts for window lengths
        n = 0:win_size(k)*fs;
        t = n*Ts;
        x = 3*sin(2*pi*0.5*t) + 2*sin(2*pi*2.5*t);
        X = abs(fft(x,N(i)));
        X = X(1:end/2);             % plot single-sided spectrum
        subplot(3,1,k);
        h = plot(f,X);
        title(strcat('{Spectrum of x(t) (win length = }', ...
            num2str(win_size(k)),')'));
        xlabel('Frequency (Hz)');
        ylabel('Spectrum');
    end
end

%% 5). FFT of Impulse Response
%
% Computing the FFT of the impulse response results in the frequency
% response of the filter. It is a low pass filter, with a cuttoff freqeuncy
% around 4.15 Hz (-3 dB).


% ----------------------- MA Filter ---------------------------
fs = 100;
M = 10;
B = ones(1,M)/M;
B1 = [B zeros(1,M)];  % increase rectangular window length
figure('Color',[1 1 1]);
subplot(2,1,1);
stem(B1);
title('h(n) = u(n)-u(n-10)');
xlabel('N');
ylabel('h(n)');
ylim([0 0.2]);

% -------------------- FFT/Frequency Response -----------------
N = 2^12;
f = linspace(0,fs/2,N/2);
X = abs(fft(B1,N));
X = X(1:end/2);
subplot(2,1,2);
h = plot(f,10*log(X));
title('FFT of h(n) (Frequency Response)');
xlabel('Frequency (Hz)');
ylabel('Spectrum (dB)');

%% 6). FFT of DC signal

DC = ones(1,fs);

N = 2^12;
f = linspace(0,fs/2,N/2);
X = abs(fft(DC,N));
X = X(1:end/2);
figure('Color',[1 1 1]);
h = plot(f,X);
title('FFT of DC signal');
xlabel('Frequency (Hz)');
ylabel('Spectrum');

%% 7). FFT of impulse signal

delta = [1 0 0 0 0 0 0 0 0 0];
N = 2^12;
f = linspace(-fs/2,fs/2,N);
X = abs(fft(delta,N));
figure('Color',[1 1 1]);
h = plot(f,X);
title('FFT of Impulse Signal');
xlabel('Frequency (Hz)');
ylabel('Spectrum');

%% 8). FFT of Cosine Signal

% --------------------- Cosine signal @ 5 Hz -----------------------
fmax = 5;
fs = 20*fmax;
Ts = 1/fs;
n = 0:1*fs;
t = n*Ts;
x = cos(2*pi*5*t);

% -------------------------- FFT ----------------------------------
N = 2^12;
f = linspace(0,fs/2,N/2);
X = abs(fft(x,N));
X = X(1:end/2);
figure('Color',[1 1 1]);
h = plot(f,X);
title('FFT of cos(2\pi5t)');
xlabel('Frequency (Hz)');
ylabel('Spectrum');

%% 9). FFT of Cosine Signal (with Gaussian Noise Added)
%
% Adding the Gaussian noise increases the peak in the spectrum at 5 Hz.
% Additionally, the rest of the spectrum trends around a magnitude of 10,
% where in the pure signal, the side lobes average ~2. 

% --------------------- Add Gaussian Noise -----------------------
n = randn(1,numel(x));
x1 = x + n;

% -------------------------- FFT ----------------------------------
N = 2^12;
f = linspace(0,fs/2,N/2);
X1 = abs(fft(x1,N));
X1 = X1(1:end/2);
figure('Color',[1 1 1]);
h = plot(f,X1);
title('FFT of cos(2\pi5t) with Added Gaussian Noise');
xlabel('Frequency (Hz)');
ylabel('Spectrum');
hold on;
h = plot(f,X,'r');
legend('cos(2\pi5t) with Added Gaussian Noise','cos(2\pi5t)');

%% 10). FFT of Signal With Three Different Cosine Components
%
% x(t) = 10cos(2pi*19t) + 3cos(2pi*2t) + 0.5cos(2pi*50t)

% ---------------------- Signal Gen -------------------------------
fmax = 50;
fs = 5*fmax;
Ts = 1/fs;
n = 0:3*fs;
t = n*Ts;
x = 10*cos(2*pi*19*t) + 3*cos(2*pi*2*t) + 0.5*cos(2*pi*50*t);
figure('Color',[1 1 1]);
subplot(2,1,1);
h = plot(t,x);
title('x(t) = 10cos(2\pi19t) + 3cos(2\pi2t) + 0.5cos(2\pi50t)');
xlabel('Time (s)');
ylabel('Amplitude');

% ----------------------- FFT ------------------------------------
N = 2^12;
X = abs(fft(x,N));
X = X(1:end/2);
f = linspace(0,fs/2,N/2);
subplot(2,1,2);
h = plot(f,X);
title('Spectrum of x(t)');
xlabel('Frequency (Hz)');
ylabel('Spectrum');

%% 11). FFT of Signal With Cosine Components Close in Frequency
%
% x(t) = 10cos(2pi*10t) + 3cos(2pi*11t) + 0.5cos(2pi*11.5t)
% 
% I used a Blackman window to remove side lobes (easier to see impulses
% close together). 

% ---------------------- Signal Gen -------------------------------
fmax = 13;
fs = 10*fmax;
Ts = 1/fs;
n = 0:10*fs;
t = n*Ts;
x = 10*cos(2*pi*10*t) + 3*cos(2*pi*11*t) + 0.5*cos(2*pi*11.5*t);
figure('Color',[1 1 1]);
subplot(2,1,1);
h = plot(t,x);
title('x(t) = 10cos(2\pi10t) + 3cos(2\pi11t) + 0.5cos(2\pi11.5t)');
xlabel('Time (s)');
ylabel('Amplitude');

% ----------------------- FFT ------------------------------------
r = 4;
x = decimate(x,r);
fs = fix(fs/r);
x = x'.*blackman(numel(x)); 
N = 2^12;
X = abs(fft(x,N));
X = X(1:end/2);
f = linspace(0,fs/2,N/2);
subplot(2,1,2);
h = plot(f,X);
title('Spectrum of x(t) - Blackman Window');
xlabel('Frequency (Hz)');
ylabel('Spectrum');


%% 12). FFT of Rectangular Signal

% ---------------------- Signal Gen -------------------------------
fs = 100;
x = [ones(1,10),zeros(1,10)];
figure('Color',[1 1 1]);
subplot(2,3,1);
h = stem(x);
title('Rectangular Signal');
ylim([0 2]);

% ----------------------- FFT ------------------------------------
N = 2^12;
X = abs(fft(x,N));
X = X(1:end/2);
f = linspace(0,fs/2,N/2);
subplot(2,3,4);
h = plot(f,X);
title({'Spectrum of'; 'Rectangular Signal'});
xlabel('Frequency (Hz)');
ylabel('Spectrum');
hold on;

%% 13). FFT of Rectangular Signal (Double Length)
%
% As signal length increases, the sidelobes get closer together.

% ---------------------- Signal Gen -------------------------------
x = ones(1,20);

subplot(2,3,2);
h = stem(x);
title({'Rectangular Signal'; '(Double Length)'});
ylim([0 2]);

% ----------------------- FFT ------------------------------------
N = 2^12;
X = abs(fft(x,N));
X = X(1:end/2);
f = linspace(0,fs/2,N/2);
subplot(2,3,5);
h = plot(f,X);
title({'Spectrum of Rectangular'; 'Signal (Double Length)'});
xlabel('Frequency (Hz)');
ylabel('Spectrum');
hold on;

%% 14). FFT of Rectangular Signal (Half Length)
%
% The opposite of the above is true here. The shorter length signal has
% sidelobes futher apart.

% ---------------------- Signal Gen -------------------------------
fs = 100;
x = [ones(1,5),zeros(1,15)];
subplot(2,3,3);
h = stem(x);
title({'Rectangular Signal'; '(Half Length)'});
ylim([0 2]);

% ----------------------- FFT ------------------------------------
N = 2^12;
X = abs(fft(x,N));
X = X(1:end/2);
f = linspace(0,fs/2,N/2);
subplot(2,3,6);
h = plot(f,X);
title({'Spectrum of Rectangular'; 'Signal (Half Length)'});
xlabel('Frequency (Hz)');
ylabel('Spectrum');

%% 15). FFT of Truncated Sinc Signal

% ---------------------- Sinc Gen -------------------------------
fs = 100
n = 0:1500;
x = sinc(t);

% --------------------- Truncate ---------------------------------
figure('Color',[1 1 1]);
subplot(2,1,1);
h = plot([x(1:600),zeros(1,900)]);
title('Truncated Sinc');
xlabel('Sample');
ylabel('Amplitude');

% ----------------------- FFT ------------------------------------
r = 4;
x = decimate(x,r);
fs = fix(fs/r);
N = 2^12;
X = abs(fft(x,N));
X = X(1:end/2);
f = linspace(0,fs/2,N/2);
subplot(2,1,2);
h = plot(f,X);
title('Spectrum of Truncated Sinc');
xlabel('Frequency (Hz)');
ylabel('Spectrum');

%% 16). FFT of Decaying Exponential Signal

% ---------------------- Exp. Decay Gen -------------------------------
fmax = 13;
fs = 10*fmax;
Ts = 1/fs;
n = 0:2*fs;
t = n*Ts;
x = exp(-1.5*t);
figure('Color',[1 1 1]);
subplot(2,1,1);
h = plot(x);
title('Decaying Exponential');
xlabel('Time (s)');
ylabel('Amplitude');

% ----------------------- FFT ------------------------------------
N = 2^12;
X = abs(fft(x,N));
X = X(1:end/2);
f = linspace(0,fs/2,N/2);
subplot(2,1,2);
h = plot(f,X);
title('Spectrum of Decaying Exponential');
xlabel('Frequency (Hz)');
ylabel('Spectrum');

%% 17). FFT of Amplitude Modulated Signal

% ---------------------- Exp. Decay Gen -------------------------------
fs = 200;
Ts = 1/fs;
n = 0:2*fs;
t = n*Ts;
x = 10*cos(2*pi*2*t).*cos(2*pi*30*t);
figure('Color',[1 1 1]);
subplot(2,1,1);
h = plot(x);
title('Amplitude Modulated Signal');
xlabel('Time (s)');
ylabel('Amplitude');

% ----------------------- FFT ------------------------------------
N = 2^12;
X = abs(fft(x,N));
X = X(1:end/2);
f = linspace(0,fs/2,N/2);
subplot(2,1,2);
h = plot(f,X);
title('Spectrum of Amplitude Modulated Signal');
xlabel('Frequency (Hz)');
ylabel('Spectrum');

%% 18). Frequency Response of FIR Filters
%
% The results of frequency response using the FFT exactly match the result
% when using freqz (from previous lab).

% -------------------- Filters Used in Previous Lab -------------------
% 10th order lowpass FIR filter with fc = 30 Hz and fs = 125 Hz;
fs = 125;
fc = 30;
Wc = fc/(fs/2);
B1 = fir1(10,Wc);

% 10th order Highpass FIR filter with fc = 30 Hz and fs = 125 Hz;
B2 = fir1(10,Wc,'high');

% 100th order Highpass FIR filter with fc = 30 Hz and fs = 125 Hz;
B3 = fir1(100,Wc,'high');

% 100th order bandpass FIR filter with fc1 = 15 Hz, fc2 = 30 Hz,
% and fs = 125 Hz;
fc1 = 15;
fc2 = 30;
Wc1 = fc1/(fs/2);
Wc2 = fc2/(fs/2);
B4 = fir1(100,[Wc1 Wc2]);

% 100th order bandstop FIR filter with fc1 = 15 Hz, fc2 = 30 Hz,
% and fs = 125 Hz;
B5 = fir1(100,[Wc1 Wc2],'stop');

BC = {B1,B2,B3,B4,B5; 
    "10th Order LP", "10th Order HP", "100th Order HP", "100th Order BP", "100th Order BS"};

% --------------------------- FFT ---------------------------------------
N = 2^12;
f = linspace(0,fs/2,N/2);
figure('Color',[1 1 1]);
sgtitle('Freqeuency Response of FIR Filters Using the FFT');

for k = 1:length(BC)
    X = abs(fft(BC{1,k},N));
    X = X(1:end/2);
    subplot(2,5,k);
    h = stem(BC{1,k});
    set(h,'MarkerSize',4);
    title(BC{2,k});
    subplot(2,5,k+5);
    h = plot(f,X);
    title({'FFT of'; BC{2,k}});
    ylim([0 1.2]);
    xlabel('Frequency (Hz)');
    ylabel('Spectrum');
end

%% 19). FFT of Real Signals
%
% * ICP
% * ECG
% * Speech
% * Audio
%
% The ICP and ECG signals have frequency components around 2 and 4 Hz.
% I believe these are heartbeat and respiration, respectivley. The ECG
% signal has a few more components up to ~22 Hz. Some of this is noise.
% 
% Both the speech and music spectrums have many freqeuncy components. The
% multitude of components shows the very high unique frequency content of
% audio. 

figure('Color',[1 1 1]);
sgtitle('Frequency Response of FIR Filters Using the FFT');

N = 2^12;

% ------------------------------ ICP -----------------------------------
load ICP.mat
fs_icp = 125;
t_icp = (0:length(icp1)-1)/fs_icp;
subplot(2,4,1);
h = plot(t_icp,icp1);
title('"ICP.mat"','FontSize',8);
xlabel('Time (s)');
ylabel('ICP (mmHg)');

% Decimate
r = 10;
icp1 = decimate(icp1,r);
fs_icp = fix(fs_icp/r);

%  Remove DC Trend
Wp = 0.1/(fs_icp/2);
[B,A] = ellip(6,0.5,20,Wp,'high');
icp1 = filtfilt(B,A,icp1);

% FFT
f_icp = linspace(0,fs_icp/2,N/2);
X_icp = abs(fft(icp1,N));
X_icp = X_icp(1:end/2);
subplot(2,4,5);
h = plot(f_icp,X_icp);
title({'FFT of'; '"ICP.mat"'},'FontSize',8);
xlabel('Frequency (Hz)');
ylabel('Spectrum');

% --------------------------------- ECG ----------------------------------
load ECGNoisy60Hz.mat
fs_ecg = 500;
t_ecg = (0:length(ecgn)-1)/fs_ecg;
subplot(2,4,2);
plot(t_ecg,ecgn);
title('"ECGNoisy60Hz.mat"','FontSize',8);
xlabel('Time (s)');
ylabel('ECG');

% Decimate
r = 10;
ecgn = decimate(ecgn,r);
fs_ecg = fix(fs_ecg/r);

%  Remove DC Trend
Wp = 0.1/(fs_ecg/2);
[B,A] = ellip(6,0.5,20,Wp,'high');
ecgn = filtfilt(B,A,ecgn);

% FFT
f_ecg = linspace(0,fs_ecg/2,N/2);
X_ecg = abs(fft(ecgn,N));
X_ecg = X_ecg(1:end/2);
subplot(2,4,6);
h = plot(f_ecg,X_ecg);
title({'FFT of'; '"ECGNoisy60Hz.mat"'},'FontSize',8);
xlabel('Frequency (Hz)');
ylabel('Spectrum');

% ---------------------------- Speech -----------------------------------
load mySpeech.mat
fs_sp = 8000;
t_sp = (0:length(speech_arr)-1)/(fs_sp);
subplot(2,4,3);
plot(t_sp,speech_arr);
title('"My Speech"','FontSize',8);
xlabel('Time (s)');
ylabel('Amplitude');

% FFT
f_sp = linspace(0,fs_sp/2,N/2);
X_sp = abs(fft(speech_arr,N));
X_sp = X_sp(1:end/2);
subplot(2,4,7);
h = plot(f_sp,X_sp);
title({'FFT of'; '"My Speech"'},'FontSize',8);
xlabel('Frequency (Hz)');
ylabel('Spectrum');


% --------------------------- Audio/Music --------------------------------
load myMusic.mat
fs_mu = 44100;
t_mu = (0:length(music_arr)-1)/(fs_mu);
subplot(2,4,4);
plot(t_mu,music_arr);
title('My Music','FontSize',8);
xlabel('Time (s)');
ylabel('Amplitude');

% FFT
f_mu = linspace(0,fs_mu/2,N/2);
X_mu = abs(fft(music_arr,N));
X_mu = X_mu(1:end/2);
subplot(2,4,8);
h = plot(f_mu,X_mu);
title({'FFT of'; '"My Music"'},'FontSize',8);
xlabel('Frequency (Hz)');
ylabel('Spectrum');
