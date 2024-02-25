function [X,f,f0] = ComputeSpectrum(x,fs,fmax,N);
% ComputeSpectrum: computes spectrum of x
%
% INPUTS
% x = input signal
% fs = sampling frequency
% fmax = max frequnecy of interest
% N = number of points to evaluate FFT
% 
% OUTPUTS
% X = spectrum of x
% f = frequencies
% f0 = frequency of highest amplitude
%
% description: 1). Resamples signal to sampling frequency 2*fmax, 
% 2). removes trend, 3). computes N-point fft
%
% Example: Compute the spectrum of an ICP signal for 0-8 Hz
% 
% load ICP.mat
% fs = 125;
% [X,f] = ComputeSpectrum(icp1,125,8,2^12);
%
%
% Ver: 1.0

% Step 1: Decimate signal to fmax
R = fix((fs/2)/fmax);
fs = fix(fs/R);
x = decimate(x,R);

% Step 2: Remove Local trend
Wp = 0.1/(fs/2);
[B,A] = ellip(6,0.5,20,Wp,'high');
x = filtfilt(B,A,x);

% Step 3: Spectrum 
X = abs(fft(x.*blackman(length(x)),N));
X = X(1:end/2);
f = linspace(0,fs/2,N/2);

figure('Color',[1 1 1]);
h = plot(f,X);
xlabel('Frequency (Hz)');
ylabel('Spectrum');

% Step 4: frequency estimation
[m,i] = max(X);
f0 = f(i);
hold on;
plot(f0,X(i),'r.');
