%% Lab 1
% by Dominic Meads

%% 1). fs = 50*fmax;

fmax = 3.5;
fs = 50*fmax;
Ts = 1/fs;
n = 0:2*fs-1;
t = n*Ts;
x = cos(2*pi*0.5*t) + cos(2*pi*1.5*t) + cos(2*pi*2.5*t) + cos(2*pi*3.5*t);

figure('Color',[1 1 1]);
h = plot(t,x);
title('x(t) Sampled at fs = 50*fmax');
ylabel('Signal');
xlabel('Time (s)');

%% 2).,3)., and 4). Plot and describe with different sampling frequencies
% each time the sampling frequency is reduced, the shape of the waveform 
% gets less accurate when compared to the signal sampled at fs = 50*fmax

figure('Color',[1 1 1]);

fs = 10*fmax;
Ts = 1/fs;
n = 0:2*fs-1;
t = n*Ts;
x = cos(2*pi*0.5*t) + cos(2*pi*1.5*t) + cos(2*pi*2.5*t) + cos(2*pi*3.5*t);
subplot(3,1,1);
h = plot(t,x);
title('x(t) Sampled at fs = 10*fmax');
ylabel('Signal');
xlabel('Time (s)');

fs = 5*fmax;
Ts = 1/fs;
n = 0:2*fs-1;
t = n*Ts;
x = cos(2*pi*0.5*t) + cos(2*pi*1.5*t) + cos(2*pi*2.5*t) + cos(2*pi*3.5*t);
subplot(3,1,2);
h = plot(t,x);
title('x(t) Sampled at fs = 5*fmax');
ylabel('Signal');
xlabel('Time (s)');

fs = 2*fmax;
Ts = 1/fs;
n = 0:2*fs-1;
t = n*Ts;
x = cos(2*pi*0.5*t) + cos(2*pi*1.5*t) + cos(2*pi*2.5*t) + cos(2*pi*3.5*t);
subplot(3,1,3);
h = plot(t,x);
title('x(t) Sampled at fs = 2*fmax');
ylabel('Signal');
xlabel('Time (s)');

%% Prefilters

figure('Color',[1 1 1]);

% no prefilter
fs = fmax;
Ts = 1/fs;
n = 0:2*fs-1;
t = n*Ts;
x = cos(2*pi*0.5*t) + cos(2*pi*1.5*t) + cos(2*pi*2.5*t) + cos(2*pi*3.5*t);
subplot(3,1,1);
h = plot(t,x);
title('x(t) Sampled at fs = fmax WITHOUT prefilter');
ylabel('Signal');
xlabel('Time (s)');
% aliased signal reconstruction (note the top two subplots are the same)
xa = 1 + cos(2*pi*0.5*t) + cos(2*pi*1*t) + cos(2*pi*1.5*t);
subplot(3,1,2);
h = plot(t,xa);
title('Aliased Signal Reconstruction');
ylabel('Signal');
xlabel('Time (s)');


% ideal prefilter
% All components with frequencies greater than 0.5*fs (1.75 Hz) are filtered out
x = cos(2*pi*0.5*t) + cos(2*pi*1.5*t);
subplot(3,1,3);
h = plot(t,x);
title('x(t) Sampled at fs = fmax WITH prefilter');
ylabel('Signal');
xlabel('Time (s)');

%% 7). Plotting ICP

load('ICP.mat');
fs = 125;
t = (0:length(icp1)-1)/(fs);
icps = icp1(2250:3000); % segment
ts = (0:length(icps)-1)/(fs);
figure('Color',[1 1 1]);
subplot(2,1,1);
plot(t,icp1);
title('Full ICP1 Signal');
xlabel('Time (s)');
ylabel('ICP (mmHg)');
subplot(2,1,2);
plot(ts,icps);
title('ICP1 Segment');
xlabel('Time (s)');
ylabel('ICP (mmHg)');

%% 8). Plotting ECG

load('ECG.mat');
fs = 500;
t = (0:length(ecg)-1)/(fs);
ecgs = ecg(20000:25000); % segment
ts = (0:length(ecgs)-1)/(fs);
figure('Color',[1 1 1]);
subplot(2,1,1);
plot(t,ecg);
title('ECG');
xlabel('Time (s)');
ylabel('Full ECG Signal');
subplot(2,1,2);
plot(ts,ecgs);
title('ECG Segment');
xlabel('Time (s)');
ylabel('ECG');

%% 9). Loading handel.mat

load('handel.mat');

fs = 8192;
t = (0:length(y)-1)/(fs);
figure('Color',[1 1 1]);
plot(t,y);
title('Handel');
xlabel('Time (s)');
ylabel('Amplitude');

%% 10). Recording my own sound

my_sound = audiorecorder(44100,16,1);
disp('Start speaking\n');
recordblocking(my_sound,2);
disp('End of recording');
p = play(my_sound);
fs = 44100;
my_sound_arr = getaudiodata(my_sound);
t = (0:length(my_sound_arr)-1)/(fs);
figure('Color',[1 1 1]);
plot(t,my_sound_arr);
title('My Sound');
xlabel('Time (s)');
ylabel('Amplitude');










