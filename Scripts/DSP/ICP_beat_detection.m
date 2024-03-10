%% project testing

% close all
clear
clc

%% loading signal

load ICP.mat
fs = 125;
t = (0:length(icp1)-1)/fs;

%% Plot Expert Annotations

% ---------------------------- ICP1 ----------------------------------

figure('Color',[1 1 1]);
sgtitle('ICP1 Annotations')
subplot(3,1,1)
h = plot(t,icp1);
hold on;
h = plot(dDT1/fs,icp1(dDT1),'r.');
set(h,'MarkerSize',18);
title('Expert DT')

subplot(3,1,2)
h = plot(t,icp1);
hold on;
h = plot(dJM1/fs,icp1(dJM1),'g.');
set(h,'MarkerSize',18);
title('Expert JM')

subplot(3,1,3)
h = plot(t,icp1);
hold on;
h = plot(dTT1/fs,icp1(dTT1),'.');
set(h,'MarkerSize',18);
set(h,'Color',[0.8 0 0.8]);
title('Expert TT')

% ---------------------------- ICP2 ----------------------------------

figure('Color',[1 1 1]);
sgtitle('ICP2 Annotations')
subplot(3,1,1)
h = plot(t,icp2);
hold on;
h = plot(dDT2/fs,icp2(dDT2),'r.');
set(h,'MarkerSize',18);
title('Expert DT')

subplot(3,1,2)
h = plot(t,icp2);
hold on;
h = plot(dJM2/fs,icp2(dJM2),'g.');
set(h,'MarkerSize',18);
title('Expert JM')

subplot(3,1,3)
h = plot(t,icp2);
hold on;
h = plot(dTT2/fs,icp2(dTT2),'.');
set(h,'MarkerSize',18);
set(h,'Color',[0.8 0 0.8]);
title('Expert TT')

%% Segment from ICP1
% include in discussion why picked a n-second window

win_time = 5;
tx = t(1:win_time*fs);

% Test segments
% x = icp1(1:win_time*fs); % 0 to 5 seconds
% x = icp1(300*fs:305*fs-1); % 300 to 305 seconds
% x = icp1(1300*fs:1305*fs-1); % 1300 to 1305 seconds  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ends with maxima
% x = icp1(2700*fs:2705*fs-1); % 2700 to 2705 seconds %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ends with maxima

% x = icp2(1:win_time*fs); % 0 to 5 seconds %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ends with maxima
% x = icp2(300*fs:305*fs-1); % 300 to 305 seconds
% x = icp2(1300*fs:1305*fs-1); % 1300 to 1305 seconds
x = icp2(2700*fs:2705*fs-1); % 2700 to 2705 seconds

figure('Color',[1 1 1]);
h = plot(tx,x);
title(strcat(num2str(win_time),' Second Segment of ICP1'));

%% Remove DC trend and high frequency noise INCREASE
%
% Min heart rate 26 bpm (0.42 Hz), fc11 = 0.41 (cuttoff as much DC, keep as
% much signal)
% Max heart rate 220 bpm (3.67 Hz), fc12 = 3.8 (keep order as small as possible) 
%
% Using MATLAB Filter design tool:
% Elliptic Bandpass filter designed using FDESIGN.BANDPASS.

% All frequency values are in Hz.
% Fs = 125;  % Sampling Frequency
% 
% Fstop1 = 0.4;     % First Stopband Frequency
% Fpass1 = 0.41;    % First Passband Frequency
% Fpass2 = 3.7;     % Second Passband Frequency
% Fstop2 = 3.8;     % Second Stopband Frequency
% Astop1 = 40;      % First Stopband Attenuation (dB)
% Apass  = 1;       % Passband Ripple (dB)
% Astop2 = 40;      % Second Stopband Attenuation (dB)

% load bp1_obj.mat % HOW TO DO IN FUNCTION

% figure('Color',[1 1 1]);
% zplane(bp1);
% figure('Color',[1 1 1]);
% freqz(bp1,2^10,fs);

% -------------------------- Filter and Plot New Signal ------------------

% x1 = filtfilt(bp1.sosMatrix,bp1.ScaleValues,x);
% 
% figure('Color',[1 1 1]);
% h = plot(tx,x);
% title(strcat(num2str(win_time),' Second Segment of ICP1'))
% hold on;
% h = plot(tx,x1+mean(x));


% NEW FILTER
% Min heart rate 26 bpm (0.42 Hz), fc11 = 0.41 (cuttoff as much DC, keep as
% much signal)
% 50/60 Hz noise fc12 = 48 Hz (keep order as small as possible) 
fc1 = 0.41;
fc2 = 48;
Wc1 = fc1/(fs/2);
Wc2 = fc2/(fs/2);
[B1,A1] = ellip(6,0.5,40,[Wc1 Wc2]);
figure('Color',[1 1 1]);
zplane(B1,A1);
figure('Color',[1 1 1]);
freqz(B1,A1,2^10,fs);

% -------------------------- Filter and Plot New Signal ------------------

x1 = filtfilt(B1,A1,x);

figure('Color',[1 1 1]);
h = plot(tx,x);
title(strcat(num2str(win_time),' Second Segment of ICP1'))
hold on;
h = plot(tx,x1+mean(x));


%%  Spectrogram of filtered signal 

r = 10;
x1d = decimate(x1,r);
fsd = fix(fs/r);

% [y1,f1,t1,p1] = spectrogram(x1d,13,[],2^12,fsd,'yaxis'); % less points faster computation
% figure('Color',[1 1 1]);
% meshc(t1,f1,p1);
% view(0,90);
% set(gca, 'YDir', 'normal');
% title('Spectrogram of  DECIMATED ICP1');
% xlabel('Time (s)');
% ylabel('Frequency (Hz)');
% colorbar
% 
% [y1,f1,t1,p1] = spectrogram(x1,125,[],2^12,fs,'yaxis'); 
% figure('Color',[1 1 1]);
% meshc(t1,f1,p1);
% view(0,90);
% set(gca, 'YDir', 'normal');
% title('Spectrogram of ICP1');
% xlabel('Time (s)');
% ylabel('Frequency (Hz)');
% colorbar

%% PSD Estimate
% Welch method (takes segments of signal -- heart rate changes over time). 

figure('Color',[1 1 1]);
[pxx,f] = pwelch(x1d,13,[],2^12,fsd);
plot(f,pxx)
xlabel('Frequency (Hz)')
ylabel('PSD (dB/Hz)')

%% Find max power frequencies 
% find maximumum power frequency, find 90% points as well. 
imax = find(pxx==max(pxx));
fheart = f(imax);
m = 0.85*max(pxx);

% split data
pxx1 = pxx(1:imax);
pxx2 = pxx(imax+1:end);

[~,i1] = min(abs(pxx1-m));
fh1 = f(i1); 

[~,i2] = min(abs(pxx2-m));
i2 = i2 + imax;
fh2 = f(i2);

%% Filter Original signal (notch filter) to only fheart 
% 
% Need this filter to always be stable so do check
%
Wc21 = fh1/(fs/2);
Wc22 = fh2/(fs/2);
[B2,A2] = ellip(2,0.5,40,[Wc21, Wc22]);
figure('Color',[1 1 1]);
zplane(B2,A2);
figure('Color',[1 1 1]);
freqz(B2,A2,2^10,fs);

% test stability for different frequencies
% fst = 125;
% ft = 0:0.01:3.8;
% bwt = 0.1/(fst/2);
% Wct = ft./(fst/2);
% tft = ones(1,numel(ft));
% 
% for k = 12:(numel(ft)-1)
%     [Bt,At] = ellip(2,0.5,40,[Wct(k)-bwt, Wct(k)+bwt]);
%     tft(k) = isstable(Bt,At);
% end
% 
% min(tft)

x2 = filtfilt(B2,A2,x);

figure('Color',[1 1 1]);
h = plot(tx,x);
title(strcat(num2str(win_time),' Second Segment of ICP1'))
set(h,'LineWidth',1.2);
set(h,'Color',[0.8 0.8 1]);
hold on;
h = plot(tx,x1+mean(x));
set(h,'LineWidth',1.2);
set(h,'Color',[0.3 0.3 1]);
hold on;
h = plot(tx,x2+mean(x));
set(h,'LineWidth',1.2);
set(h,'Color',[1 0.6 0.6]);
hold on;

%% Detection 

x3 = [diff(x2)', 0];
h = plot(tx,x3+mean(x));
set(h,'LineWidth',1.2);
set(h,'Color',[0.2 0.73 1]);
hold on;

a = x3;
L = length(x3);
idmax = find((a(1:L-2)<a(2:L-1))&(a(2:L-1)>a(3:L)))+1;
idmin = find((a(1:L-2)>a(2:L-1))&(a(2:L-1)<a(3:L)))+1;

h = plot(idmax/fs,a(idmax)+mean(x),'r.');
set(h,'MarkerSize',18);
hold on;
h = plot(idmin/fs,a(idmin)+mean(x),'b.');
set(h,'MarkerSize',18);
% legend('x','x1 (fundamental)','x2 (fc = max PSD)','diff(x2)','idmax','idmin');
hold on;

% intrv_srt = idmax(1);
% intrv_end = idmin(1+1);
% x_intrvl = x(intrv_srt:intrv_end);
% p1 = find(x_intrvl==max(x_intrvl),1);
% p[k] = p1+intrv_srt;

p = zeros(1,numel(idmax));

for k = 1:numel(idmax)
    intrv_srt = idmax(k);
    intrv_end = idmin(k+1);
    x_intrvl = x(intrv_srt:intrv_end);
    p1 = find(x_intrvl==max(x_intrvl),1);
    p(k) = p1+intrv_srt;
end

h = plot(p/fs,x(p),'g.');
set(h,'MarkerSize',18);





 