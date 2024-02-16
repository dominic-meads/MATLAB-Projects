%% EE430 Lab 4
% by Dominic Meads

close all;
clear;
clc;
%% Part A: Designing and Analyzing Filters

%% Lowpass FIR
% When comparing Frequency responses, the 100th order filter has a much
% steeper roll-off. 

% 10th order lowpass FIR filter with fc = 30 Hz and fs = 125 Hz;
fs = 125;
fc = 30;
Wc = fc/(fs/2);
B = fir1(10,Wc);

% Filter Test
Ts = 1/fs;
n = 0:0.5*fs; % plot for 1/2 sec
t = n*Ts;
% Test signal with frequency components @ 5 Hz and 40 Hz
x = 5*cos(2*pi*5*t) + 0.5*cos(2*pi*40*t); 
figure('Color',[1 1 1]);
h = plot(t,x);
title('Test of 10th Order LP FIR With fc = 30 Hz');
xlabel('Time (s)');
ylabel('Signal');
set(h,'LineWidth',2);
hold on;
y = filter(B,1,x);
h = plot(t,y);
legend('Original Signal','Filtered Signal');
set(h,'Color',[1 0.5 0.5]);
set(h,'LineWidth',2);
axis tight;

figure('Color',[1 1 1]);
stem(B);
title('10th Order LP Filter Coefficients');
axis tight;

figure('Color',[1 1 1]);
freqz(B,1,2^10,fs);
title('10th Order LP Frequency Response');



% 100th order lowpass FIR filter with fc = 30 Hz and fs = 125 Hz;
B = fir1(100,Wc);

figure('Color',[1 1 1]);
stem(B);
title('100th Order LP Filter Coefficients');
axis tight;


figure('Color',[1 1 1]);
freqz(B,1,2^10,fs);
title('100th Order LP Frequency Response');

%% Highpass FIR
% Similar to the lowpass FIR filter, a higher order filter results in a
% steeper (better) roll-off;    

% 10th order Highpass FIR filter with fc = 30 Hz and fs = 125 Hz;
B = fir1(10,Wc,'high');

% Filter Test
figure('Color',[1 1 1]);
h = plot(t,x);
title('Test of 10th Order HP FIR With fc = 30 Hz');
xlabel('Time (s)');
ylabel('Signal');
set(h,'LineWidth',2);
hold on;
y = filter(B,1,x);
h = plot(t,y);
legend('Original Signal','Filtered Signal');
set(h,'Color',[1 0.5 0.5]);
set(h,'LineWidth',2);
axis tight;

figure('Color',[1 1 1]);
stem(B);
title('10th Order HP Filter Coefficients');
axis tight;

figure('Color',[1 1 1]);
freqz(B,1,2^10,fs);
title('10th Order HP Frequency Response');



% 100th order highpass FIR filter with fc = 30 Hz and fs = 125 Hz;
B = fir1(100,Wc,'high');

figure('Color',[1 1 1]);
stem(B);
title('100th Order HP Filter Coefficients');
axis tight;


figure('Color',[1 1 1]);
freqz(B,1,2^10,fs);
title('100th Order HP Frequency Response');

%% Bandpass and Bandstop FIR Using "fir1"

%--------------------Bandpass---------------------------------------

% 100th order bandpass FIR filter with fc1 = 15 Hz, fc2 = 30 Hz,
% and fs = 125 Hz;
fc1 = 15;
fc2 = 30;
Wc1 = fc1/(fs/2);
Wc2 = fc2/(fs/2);
B = fir1(100,[Wc1 Wc2]);

% Filter Test
n = 0:1*fs; % plot for 1 sec
t = n*Ts;
% Test signal with frequency components @ 5 Hz, 22.5 Hz and 40 Hz
x = 10*cos(2*pi*5*t) + 4*cos(2*pi*22.5*t) + 1*cos(2*pi*40*t); 
figure('Color',[1 1 1]);
h = plot(t,x);
title('Test of 100th Order BP FIR With fc1 = 15 Hz, fc2 = 30 Hz');
xlabel('Time (s)');
ylabel('Signal');
set(h,'LineWidth',2);
hold on;
y = filter(B,1,x);
h = plot(t,y);
legend('Original Signal','Filtered Signal');
set(h,'Color',[1 0.5 0.5]);
set(h,'LineWidth',2);
axis tight;

figure('Color',[1 1 1]);
stem(B);
title('100th Order BP Filter Coefficients');
axis tight;

figure('Color',[1 1 1]);
freqz(B,1,2^10,fs);
title('100th Order BP Frequency Response');


%--------------------Bandstop---------------------------------------

% 100th order bandstop FIR filter with fc1 = 15 Hz, fc2 = 30 Hz,
% and fs = 125 Hz;
B = fir1(100,[Wc1 Wc2],'stop');

figure('Color',[1 1 1]);
h = plot(t,x);
title('Test of 100th Order BS FIR With fc1 = 15 Hz, fc2 = 30 Hz');
xlabel('Time (s)');
ylabel('Signal');
set(h,'LineWidth',2);
hold on;
y = filter(B,1,x);
h = plot(t,y);
legend('Original Signal','Filtered Signal');
set(h,'Color',[1 0.5 0.5]);
set(h,'LineWidth',2);
axis tight;

figure('Color',[1 1 1]);
stem(B);
title('100th Order BS Filter Coefficients');
axis tight;

figure('Color',[1 1 1]);
freqz(B,1,2^10,fs);
title('100th Order BS Frequency Response');

%% Bandpass FIR Using "fir2"

% 100th order bandpass FIR filter with fc1 = 15 Hz, fc2 = 30 Hz,
% and fs = 125 Hz;

fc1 = 15;
fc2 = 30;
Wc1 = fc1/(fs/2);
Wc2 = fc2/(fs/2);
F = [0 Wc1 Wc2 1]; % Frequency breakpoints (must begin with 0 and end with 1)
M = [0 1 1 0];     % Magnitude breakpoints
B = fir2(100,F,M);

% Filter Test
figure('Color',[1 1 1]);
h = plot(t,x);
title('Test of 100th Order BP FIR (fir2) With fc1 = 15 Hz, fc2 = 30 Hz');
xlabel('Time (s)');
ylabel('Signal');
set(h,'LineWidth',2);
hold on;
y = filter(B,1,x);
h = plot(t,y);
legend('Original Signal','Filtered Signal');
set(h,'Color',[1 0.5 0.5]);
set(h,'LineWidth',2);
axis tight;

figure('Color',[1 1 1]);
stem(B);
title('100th Order BP (fir2) Filter Coefficients');
axis tight;

figure('Color',[1 1 1]);
freqz(B,1,2^10,fs);
title('100th Order BP (fir2) Frequency Response');

%% Part B: Biomedical Signal Processing

load ICP.mat;
fs = 125;
tfull = (0:length(icp1)-1)/fs;
figure('Color',[1 1 1]);
h = plot(tfull,icp1);
title('ICP.mat');
xlabel('Time (s)');
ylabel('ICP (mmHg)');

icps = icp1(10^4:10^4+300); % small segment
t = (0:length(icps)-1)/fs;

%% Lowpass Filter to Smooth Quantization Noise
% I decided to cut all frequencies above 15 Hz. In one of the lab
% tutorials, we used a moving average filter to do the smoothing. By
% plotting the frequency response of the moving average filter, I found the
% cuttoff frequency equal to 15 Hz. The moving average filter resulted in a
% nice smoothing of the ICP signal, so I tried it with an FIR filter. This 
% cut out the quantization noise, while still retaining most of the
% information of the signal. 
%
% A higher order filter (n = 100) did not do much more smoothing than 30th
% order filter, but any orders below 30 and the signal starts to attenuate.
% I chose n = 30 to decrease computational load of my filter while still
% providing similar results to a 100th order filter.

fc = 15;
Wc = fc/(fs/2);
B = fir1(30,Wc);
icp_filt = filtfilt(B,1,icp1);
icp_filts = icp_filt(10^4:10^4+300); % segment of filtered signal

figure('Color',[1 1 1]);
h = plot(t,icps);
title('Segment of ICP Signal');
xlabel('Time (s)');
ylabel('ICP (mmHg)');
set(h,'LineWidth',1.5);
hold on;
h = plot(t,icp_filts);
legend('ICP signal', 'Smoothed ICP signal');
set(h,'Color',[1 0.5 0.5]);
set(h,'LineWidth',1.5);

%% FIR Highpass to Eliminate Low-frequency Trends
% The intracranial pressure of this patient cycles at about 2.5 Hz. The DC
% component is at 0 Hz. These frequencies are very close together, so I
% chose a cuttoff at 1 Hz, with a very steep roll-off (high order filter)
% to eliminate only the DC trend. 

fc = 1;
Wc = fc/(fs/2);
B = fir1(120,Wc,'high');
icp_filt = filtfilt(B,1,icp1);
icp_filts = icp_filt(10^4:10^4+300); % segment of filtered signal

figure('Color',[1 1 1]);
h = plot(t,icps);
title('Segment of ICP Signal');
xlabel('Time (s)');
ylabel('ICP (mmHg)');
set(h,'LineWidth',1.5);
ylim([-1 10]);
hold on;
h = plot(t,icp_filts);
legend('ICP signal', 'ICP signal w/o DC Component');
set(h,'Color',[1 0.5 0.5]);
set(h,'LineWidth',1.5);

%% FIR Bandpass to Eliminate DC Value and Quantization Noise
% I used the previous cuttoff frequencies (15 hz for quantization noise and
% 1 Hz for DC value) to create my bandpass filter. I kept the order high
% (n = 120) to ensure a steep roll-off. This eliminated the DC portion
% while passing the 2.5 Hz ICP fundamental frequency (close together in 
% terms of frequency). 

fc1 = 1;
fc2 = 15;
Wc1 = fc1/(fs/2);
Wc2 = fc2/(fs/2);
B = fir1(120,[Wc1 Wc2]);

icp_filt = filtfilt(B,1,icp1);
icp_filts = icp_filt(10^4:10^4+300); % segment of filtered signal

figure('Color',[1 1 1]);
h = plot(t,icps);
title('Segment of ICP Signal');
xlabel('Time (s)');
ylabel('ICP (mmHg)');
set(h,'LineWidth',1.5);
ylim([-1 10]);
hold on;
h = plot(t,icp_filts);
legend('ICP signal', 'ICP signal (No DC Value or Quant. Noise)');
set(h,'Color',[1 0.5 0.5]);
set(h,'LineWidth',1.5);

%% FIR Bandpass for Cardiac Component
% I kept fc1 = 1 Hz to eliminate the DC component (found previously). fc2
% was set to 4 Hz. This gave good balance between sinusoidal shape and
% amplitude, meaning the resulting waveform is very sinusoidal without to
% much attenuation. An increase in fc2 deforms the signal, while a decrease
% causes the amplitude to be smaller. 
%
% Additionally, I had to increase the filter order for 120 to 150 in order
% to keep the DC component at zero. I am not sure why this is, but maybe
% the fundamental frequency has a low-frequency trend very close to 1 Hz?
%
% Answer ^: Yes, increasing fc1 to 1.4 allows the order of the filter to
% drop back down to 120 (previous filters), meaning the cardiac compnent
% does have a low-frequency trend close to 1 that did not show up in the
% previous bandpass filter example. 

fc1 = 1;
fc2 = 4;
Wc1 = fc1/(fs/2);
Wc2 = fc2/(fs/2);
B = fir1(150,[Wc1 Wc2]);

icp_filt = filtfilt(B,1,icp1);
icp_filts = icp_filt(10^4:10^4+300); % segment of filtered signal

figure('Color',[1 1 1]);
h = plot(t,icps);
title('Segment of ICP Signal');
xlabel('Time (s)');
ylabel('ICP (mmHg)');
set(h,'LineWidth',1.5);
ylim([-1 10]);
hold on;
h = plot(t,icp_filts);
legend('ICP signal', 'Cardiac Component');
set(h,'Color',[1 0.5 0.5]);
set(h,'LineWidth',1.5);

%% Eliminate Noise from "ECGNoisy60Hz"
% The ECG signal is filtered at fc = 45 hz with an order of 30 to get rid
% of the noise. Peak amplitude is reduced from original signal, however
% noise is no longer present. Peak amplitude reduction can be avoided by 
% using a very high order filter (n = 300) and a lower fc.

load ECGNoisy60Hz.mat
ecgs = ecgn(10^4:10^4+750); % small segment
t = (0:length(ecgs)-1)/fs;

figure('Color',[1 1 1]);
tfull = (0:length(ecgn)-1)/fs;
plot(tfull,ecgn);
title('ECGNoisy60Hz.mat');
xlabel('Time (s)');
ylabel('ECG');

fc = 45;
Wc = fc/(fs/2);
B = fir1(200,Wc);
ecg_filt = filtfilt(B,1,ecgn);
ecg_filts = ecg_filt(10^4:10^4+750); % segment of filtered signal

figure('Color',[1 1 1]);
h = plot(t,ecgs);
title('Segment of ECG Signal');
xlabel('Time (s)');
ylabel('ECG');
set(h,'LineWidth',1.5);
hold on;
h = plot(t,ecg_filts);
legend('ECG signal', 'ECG signal (No 60 Hz Noise)');
set(h,'Color',[1 0.5 0.5]);
set(h,'LineWidth',1.5);

%% Eliminate Quantization Noise from "ECGQuantization"

load ECGQuantization1.mat
ecgQs = ECGQuantization(10^4:10^4+300); % small segment
t = (0:length(ecgQs)-1)/fs;

figure('Color',[1 1 1]);
tfull = (0:length(ECGQuantization)-1)/fs;
plot(tfull,ECGQuantization);
title('ECGQuantization1.mat');
xlabel('Time (s)');
ylabel('ECG');

fc = 23;
Wc = fc/(fs/2);
B = fir1(120,Wc);
ecgQ_filt = filtfilt(B,1,ECGQuantization);
ecgQ_filts = ecgQ_filt(10^4:10^4+300); % segment of filtered signal

figure('Color',[1 1 1]);
h = plot(t,ecgQs);
title('Segment of ECGQuantization Signal');
xlabel('Time (s)');
ylabel('ECG');
set(h,'LineWidth',1.5);
hold on;
h = plot(t,ecgQ_filts);
legend('ECG signal', 'ECG signal (No Quantization Noise)');
set(h,'Color',[1 0.5 0.5]);
set(h,'LineWidth',1.5);

%% Eliminate Baseline Drift Noise Due to Patient Movement

load ECGBaselineDrift.mat
ecgBDs = ecgn(10^4:10^4+5000); % small segment
t = (0:length(ecgBDs)-1)/fs;

figure('Color',[1 1 1]);
tfull = (0:length(ecgn)-1)/fs;
plot(tfull,ecgn);
title('ECGBaselineDrift.mat');
xlabel('Time (s)');
ylabel('ECG');

fc = 2;
Wc = fc/(fs/2);
B = fir1(150,Wc,'high');
ecgBD_filt = filtfilt(B,1,ecgn);
ecgBD_filts = ecgBD_filt(10^4:10^4+5000); % segment of filtered signal

figure('Color',[1 1 1]);
h = plot(t,ecgBDs);
title('Segment of ECG Signal');
xlabel('Time (s)');
ylabel('ECG');
ylim([0 14]);
set(h,'LineWidth',1.5);
hold on;
h = plot(t,ecgBD_filts);
legend('ECG signal', 'ECG signal (No Baseline Drift Noise)');
set(h,'Color',[1 0.5 0.5]);
set(h,'LineWidth',1.5);

%% Eliminate All Noise (60 Hz, Quantization, Baseline Drift)

load ECGBaselineCombined.mat
ecgBC = ecgnc(570000:580000); % small segment
t = (0:length(ecgBC)-1)/fs;

figure('Color',[1 1 1]);
tfull = (0:length(ecgnc)-1)/fs;
plot(tfull,ecgnc);
title('ECGBaselineCombined.mat (Not Filtered)');
xlabel('Time (s)');
ylabel('ECG');

% ---------------------- Remove Baseline Drift ----------------------
fc1 = 3;
Wc1 = fc1/(fs/2);
B1 = fir1(200,Wc1,'high');
ecgBD_filt = filtfilt(B1,1,ecgnc);
ecgBD_filts = ecgBD_filt(570000:580000); % segment of filtered signal

% figure;
% freqz(B1,1,2^10,fs);

% plot a segment
figure('Color',[1 1 1]);
h = plot(t,ecgBC);
title('Baseline Drift Removal in ECG');
xlabel('Time (s)');
ylabel('ECG');
ylim([-2 14]);
set(h,'LineWidth',0.7);
hold on;
h = plot(t,ecgBD_filts);
legend('ECG signal', 'ECG signal (No Baseline Drift)');
set(h,'Color',[0.3 0.9 0.3]);
set(h,'LineWidth',0.7);

% ---------------------- Remove 60 Hz Noise -----------------------
fc2 = 55;
Wc2 = fc2/(fs/2);
B2 = fir1(200,Wc2); % Steep roll off
ecgBD_60_filt = filtfilt(B2,1,ecgBD_filt);
ecgBD_60_filts = ecgBD_60_filt(1198750:1201250); % segment of filtered signal
t2 = (0:length(ecgBD_60_filts)-1)/fs;

% figure;
% freqz(B2,1,2^10,fs);

% plot a segment
figure('Color',[1 1 1]);
h = plot(t2,ecgnc(1198750:1201250));
title('60 Hz Noise Removal in ECG');
xlabel('Time (s)');
ylabel('ECG');
ylim([-2 14]);
set(h,'LineWidth',0.7);
hold on;
h = plot(t2,ecgBD_filt(1198750:1201250));
set(h,'Color',[0.3 0.9 0.3]);
set(h,'LineWidth',0.7);
hold on;
h = plot(t2,ecgBD_60_filts);
set(h,'Color',[1 0.5 0.5]);
set(h,'LineWidth',0.7);
legend('ECG signal', 'ECG signal (No Baseline Drift)', 'ECG signal (No 60 Hz Noise)');

% ---------------------- Remove Quantization Noise -----------------------

% Moving average filter
M = 10;
B3 = 1/M*ones(M,1); 

ecgBD_60_Q_filt = filtfilt(B3,1,ecgBD_60_filt);
ecgBD_60_Q_filts = ecgBD_60_Q_filt(1200075:1200400); % segment of filtered signal
t3 = (0:length(ecgBD_60_Q_filts)-1)/fs;

% figure;
% freqz(B3,1,2^10,fs);

% plot a segment
figure('Color',[1 1 1]);
h = plot(t3,ecgBD_60_filt(1200075:1200400));
title('Quantization Noise Removal in ECG');
xlabel('Time (s)');
ylabel('ECG');
ylim([-0.6 2.1]);
set(h,'LineWidth',1.5);
set(h,'Color',[1 0.5 0.5]);
hold on;
h = plot(t3,ecgBD_60_Q_filts);
set(h,'Color',[0.8 0.3 0.8]);
set(h,'LineWidth',1.5);
legend('ECG signal (No Baseline Drift/60 Hz Noise)','ECG Signal (No Quantization Noise');

% ---------------------- Plot the Filtered Signal ----------------------

% plot entire signal
figure('Color',[1 1 1]);
h = plot(tfull,ecgnc);
title('ECG Signal With and Without Noise');
xlabel('Time (s)');
ylabel('ECG');
ylim([-2 14]);
set(h,'LineWidth',0.7);
hold on;
h = plot(tfull,ecgBD_60_Q_filt);
set(h,'Color',[0.8 0.3 0.8]);
set(h,'LineWidth',0.7);
legend('ECG signal', 'ECG Signal (Filtered)');

% plot a segment
ecg_orig_seg = ecgnc(1199750:1201500);
ecg_filtered_seg = ecgBD_60_Q_filt(1199750:1201500);
t4 = (0:length(ecg_filtered_seg )-1)/fs;
figure('Color',[1 1 1]);
h = plot(t4,ecg_orig_seg);
title('ECG Segments With and Without Noise');
xlabel('Time (s)');
ylabel('ECG');
ylim([-0.8 10]);
set(h,'LineWidth',0.7);
hold on;
h = plot(t4,ecg_filtered_seg);
legend('Original ECG signal', 'ECG signal (Filtered)');
set(h,'Color',[0.8 0.3 0.8]);
set(h,'LineWidth',0.7);

%% Part C: Speech Processing

%% Recording and Plotting My Speech
% When I speak, a microphone sensor and analog circuit turns my speech into 
% a continuous analog signal. In order to avoid aliasing, the analog signal
% is pass through a lowpass filter with a cuttoff frequency equal to half
% of the sampling frequency. The analog filtered signal is then sampled and 
% held over a certain interval (sampling period) which turns the continous 
% signal into a discrete time signal. This is done by a sample/hold circuit
% or the ADC. The discrete time signal is then quantized (rounded to 
% a certain number of bits and turned into a digital number) by the ADC. 
% This digital quantity is the one we see plotted. 

load mySpeech.mat
fs = 8000;
sound(speech_arr,fs);
pause(10);
% speech = audiorecorder(fs,16,1);
% disp('Start speaking\n');
% recordblocking(speech,10);
% disp('End of recording');
% p = play(speech);
% speech_arr = getaudiodata(speech);
% t = (0:length(speech_arr)-1)/(fs);
% save("mySpeech.mat","speech_arr","t");
figure('Color',[1 1 1]);
plot(t,speech_arr);
title('My Speech');
xlabel('Time (s)');
ylabel('Amplitude');

%% Applying Lowpass Filters to Speech
% As the cuttoff frequency decreases, the speech sample get more muffled
% and sounds lower pitched. This makes sense, because the low pass filters
% are removing some of the high pitch (high frequency) components of the
% speech. 

fc = [3000 2500 1500 1000];
Wc = zeros(4,1);
order = 30;
B = zeros(4,(order+1));
y = zeros(4,length(speech_arr));

for k = 1:length(fc)
    Wc(k) = fc(k)/(fs/2);
    B(k,:) = fir1(order,Wc(k));
    y(k,:) = filtfilt(B(k,:),1,speech_arr);
    sound(y(k,:),fs);
    pause(10);
end

%% Applying Highpass Filters to Speech
% The recordings increase in pitch as the cuttoff frequency increases,
% because the low frequency components are being attenuated by the filter.
% When fc = 2000 Hz, I cannot hear the recording well, meaning it does not
% have many components above 2000 Hz. This makes sense because I have a
% lower voice. 

fc = [100 500 1000 2000];
Wc = zeros(4,1);
order = 30;
B = zeros(4,(order+1));
y = zeros(4,length(speech_arr));

for k = 1:length(fc)
    Wc(k) = fc(k)/(fs/2);
    B(k,:) = fir1(order,Wc(k),'high');
    B1 = B(k,:);
    y(k,:) = filtfilt(B(k,:),1,speech_arr);
    sound(y(k,:),fs);
    pause(10);
end

%% Adding Gaussian Noise
% I tried a moving average filter to smooth the signal noise (similar to 
% eliminating quantization noise). The signals sounded best with M = 10,
% and the words can be made out in the first signal (where noise is at 0.01
% amplitude). In the higher noise signals, I can hear speech, but the words
% are not discernable.

% --------------------------- Add Noise ------------------------------
r = randn(3,length(speech_arr));
ygauss = zeros(3,length(speech_arr));
figure('Color',[1 1 1]);
subplot(4,1,1);
h = plot(t,speech_arr);
title('Speech');
ylabel('Amplitude');
xlabel('Time (s)');
for k = 1:3
    ygauss(k,:) = speech_arr' + (0.01*k).*r(k,:); % Add different amounts of noise
    % sound(ygauss(k,:),fs);
    % pause(10);
    subplot(4,1,k+1)
    h = plot(t,ygauss(k,:));
    title(strcat('{Speech Signal With }', num2str(0.01*k), '*randn Added'));
    ylabel('Amplitude');
    xlabel('Time (s)');
end 


% -------------------------- Try to Filter ---------------------------- 

M = 10;
B1 = 1/M*ones(M,1);

yfilt = zeros(3,length(speech_arr));
figure('Color',[1 1 1]);
subplot(4,1,1);
h = plot(t,speech_arr);
title('Speech');
ylabel('Amplitude');
xlabel('Time (s)');
for k = 1:3
    yfilt(k,:) = filtfilt(B1,1,ygauss(k,:)); 
    sound(yfilt(k,:),fs);
    pause(10);
    subplot(4,1,k+1)
    h = plot(t,yfilt(k,:));
    title(strcat('{Speech Signal With }', num2str(0.01*k), ['*randn Added' ...
        ' (Then Passed Through Moving Average Filter)']));
    ylabel('Amplitude');
    xlabel('Time (s)');
end 

%% Adding Bandlimited Gaussian Noise
% Choose Bandwidth = 1.25 KHz
% Noise1 is 1-1.25KHz
% Noise2 is 1.25-2.5KHz
% Noise3 is 2.5-3.75KHz
%
% I tried to use stopband filters with a high order (n = 200) to filter out
% the noise in my signal. Each stopband filter cuttoff frequency
% corresponds to the cuttoff frequency of bandlimited Gaussian noise. The
% results are okay. I can hear that there is speech behind the noise,
% however is is difficult to make out. 

% -------------------------- Adding Noise ------------------------
r = randn(3,length(speech_arr));
rband = zeros(3,length(speech_arr));
fc = [1 1250 2500 3750];
Wc = zeros(4,1);
order = 30;
B = zeros(3,(order+1));
y = zeros(3,length(speech_arr));
figure('Color',[1 1 1]);
subplot(4,1,1);
h = plot(t,speech_arr);
title('Speech');
ylabel('Amplitude');
xlabel('Time (s)');
for k = 1:3
    Wc(k) = fc(k)/(fs/2);
    Wc(k+1) = fc(k+1)/(fs/2);
    B(k,:) = fir1(order,[Wc(k) Wc(k+1)]);
    rband(k,:) = filtfilt(B(k,:),1,r(k,:)); % bandlimit randn noise
    y(k,:) = rband(k,:) + speech_arr';      % add banded noise 
    sound(y(k,:),fs);
    pause(10);
    subplot(4,1,k+1)
    h = plot(t,y(k,:));
    title(strcat('{Speech Signal With Bandlimited Noise Added (fc1 = }', num2str(fc(k)), '{ and fc2 = }', num2str(fc(k+1)), ')')); 
    ylabel('Amplitude');
    xlabel('Time (s)');
end

% -------------------------- Filtering Noise ------------------------

order = 200;
B = zeros(3,(order+1));
yfilt = zeros(3,length(speech_arr));
figure('Color',[1 1 1]);
subplot(4,1,1);
h = plot(t,speech_arr);
title('Speech');
ylabel('Amplitude');
xlabel('Time (s)');
for k = 1:3
    Wc(k) = fc(k)/(fs/2);
    Wc(k+1) = fc(k+1)/(fs/2);
    B(k,:) = fir1(order,[Wc(k) Wc(k+1)],'stop');
    yfilt(k,:) = filtfilt(B(k,:),1,y(k,:));   
    sound(yfilt(k,:),fs);
    pause(10);
    subplot(4,1,k+1)
    h = plot(t,yfilt(k,:));
    title(strcat('{Noisy Speech Signal Filtered with Bandstop (fc1 = }', num2str(fc(k)), '{ and fc2 = }', num2str(fc(k+1)), ')')); 
    ylabel('Amplitude');
    xlabel('Time (s)');
end


%% Adding Sinusoidal Interference
% I was able to use a high order notch filter to eliminate some of the
% noise. The speech is easily understandable, however I needed a very high
% order filter (n = 700) to attenuate the stopband @ 3000 Hz. 

% ------------------------ Adding Noise ---------------------------
fs = 8000;
Ts = 1/fs;
n = 0:10*fs-1;             % 10 seconds
t = n*Ts;
x = 0.5*cos(2*pi*3000*t);  % 3000 Hz noise
y = speech_arr'+x;         % add noise
figure('Color',[1 1 1]);
subplot(3,1,1);
h = plot(t,speech_arr);
title('Speech');
ylabel('Amplitude');
xlabel('Time (s)');
subplot(3,1,2);
h = plot(t,y);
title('Speech With 3000 Hz Sinusoidal Noise');
ylabel('Amplitude');
xlabel('Time (s)');
sound(y,fs);
pause(10);

% ------------------------ Filtering ---------------------------

fnotch = 3000;                      % notch filter cutoff
band = 10;                          % Stopband width
fc = [fnotch-band, fnotch+band];    % create stopband
Wc = fc/(fs/2);
B = fir1(700,Wc,'stop');            % high order for stopband attenuation
yfilt = filtfilt(B,1,y);

subplot(3,1,3);
h = plot(t,yfilt);
title('Noisy Speech Filtered with Bandstop FIR Filter (fc1 = 2990, fc2 = 3010)');
ylabel('Amplitude');
xlabel('Time (s)');
sound(yfilt,fs);
pause(10);

figure('Color',[1 1 1]);
freqz(B,1,2^10,fs);

%% Part D: Audio Processing and Digital Equalization

%% Loading Music
% The digital sound signal in bit form (just a number) is sent over certain
% intervals (sampling period) to a DAC (digital to analog converter). The
% DAC takes the digital number, and converts it to an analog voltage. The
% analog voltage signal is then passed though a lowpass filter to get rid
% of any images in its spectrum (anti-imaging filter). This helps prevent
% unwanted aliasing of the final signal. Additionally, the lowpass filter
% can help smooth the signal. Once past the low pass filter, the smoothed
% analog signal is pushed through an audio amplifier and speaker, which
% allows us to hear the music. 

fs = 44100;
load myMusic.mat
sound(music_arr,fs);
pause(10);
% music = audiorecorder(fs,16,1);
% disp('Start playing music\n');
% recordblocking(music,10);
% disp('End of recording');
% p = play(music);
% music_arr = getaudiodata(music);
% t = (0:length(music_arr)-1)/(fs);
% save("myMusic.mat","music_arr","t");
figure('Color',[1 1 1]);
plot(t,music_arr);
title('My Music');
xlabel('Time (s)');
ylabel('Amplitude');

%% Designing an Equalizer with "fir1"
% There is some discrepancy in the frequency ranges of bass, mid, and
% treble, so mine are defined as:
% Bass signals are from 20 - 200 Hz,
% Mid signals are from 200 - 5K Hz,
% Treble signals are from 5k - 20K Hz.

fc = [20 200 5000 20000];
Wc = fc/(fs/2);

bass = fir1(100,[Wc(1) Wc(2)]);
figure('Color',[1 1 1]);
freqz(bass,1,2^10,fs);
title('Bass Bandpass Filter');
mid  = fir1(100,[Wc(2) Wc(3)]);
figure('Color',[1 1 1]);
freqz(mid,1,2^10,fs);
title('Mid Bandpass Filter');
treble = fir1(100,[Wc(3) Wc(4)]);
figure('Color',[1 1 1]);
freqz(treble,1,2^10,fs);
title('Treble Bandpass Filter');

%% Designing an Equalizer with "fir2"
% Input frequency ranges and magnitudes for three bands.

b_low = 0;     % must be 0
b_high = 200;
b_mag = 0;

m_low = 201;
m_high = 5000;
m_mag = 2;

t_low = 5001;
t_high = fs/2; % must be fs/2
t_mag = 0.5;

f = [b_low, b_high, m_low, m_high, t_low, t_high]/(fs/2);
m = [b_mag, b_mag, m_mag, m_mag, t_mag, t_mag];
b = fir2(800,f,m);
figure('Color',[1 1 1]);
freqz(b,1,2^10,fs);
title('Equalizer Multi-Band Filter');

%% Applying the filters
% Both filter types designed using "fir1" and "fir2" work, however the
% magnitude adjustment of the multi-band filter using "fir2" is more fun to
% use. I can control how much of what kind of sound gets though, and play
% with the music. With the "fir1" filters, they are just bandpass, so the
% frequencies either get through, or they dont. 

% Use fir1 bandpass filter to filter out mid and treble;
music_filt1 = filtfilt(bass,1,music_arr);
figure('Color',[1 1 1]);
plot(t,music_filt1);
title('My Music With No Mid or Treble');
xlabel('Time (s)');
ylabel('Amplitude');
sound(music_filt1,fs);
pause(10);

% Use fir2 bandpass filter to eliminate bass, attenuate some mid, and amplify treble;
m = [0 0 0.3 0.3 2 2]; % set magnitudes
b = fir2(800,f,m);
music_filt2 = filtfilt(b,1,music_arr);
figure('Color',[1 1 1]);
plot(t,music_filt2);
title('My Music With No Bass, Attenuated Mid, and Amplified Treble');
xlabel('Time (s)');
ylabel('Amplitude');
sound(music_filt2,fs);
pause(10);



















