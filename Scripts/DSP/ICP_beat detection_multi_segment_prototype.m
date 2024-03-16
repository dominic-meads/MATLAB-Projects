%% ICP Function prototyping 

close all
clear
clc

% --------- initialize matrix to hold 10 sec segments of ICP signal -----
load ICP.mat
% sig = icp2(150*125+1:180*125); % 40 sec test segment
% sig = icp1(1:40*125); % 20 sec test segment
% sig = icp1(1760*125+1:1810*125); 
sig = icp1;
% sig = icp2;
len = length(sig);
fs = 125;
t = (0:length(sig)-1)/fs;
win_time = 10;
seg_len = win_time*fs;
num_seg = ceil(length(sig)/seg_len);
icp_seg = zeros(num_seg,seg_len);

% --------------------- fill matrix with icp signal --------------------

for k = 1:num_seg
    i_start = (k-1)*seg_len+1;               % start index
    i_stop = k*seg_len;                      % stop index
    if i_stop <= len                         % fill matrix with the entire segment
        icp_seg(k,:) = sig(i_start:i_stop);
    else                                     % fill matrix with remainder of signal
        rem = len-((k-1)*seg_len);
        icp_seg(k,1:rem) = sig(i_start:len);
    end 
end

% ----------------- Filter low freq trend and high freq noise ---------- 
% Min heart rate 27 bpm (0.45 Hz), fc1_1 = 0.45 
% link: 
% Max heart rate (Ventricular Tachycardia) 250 bpm (4.17 Hz), fc1_2 = 4.17
% link: https://ecgwaves.com/topic/ventricular-tachycardia-vt-ecg-treatment-causes-management/

fc1_1 = 0.45;
fc1_2 = 4.17;
Wc1_1 = fc1_1/(fs/2);
Wc1_2 = fc1_2/(fs/2);
[B1,A1] = ellip(4,0.5,40,[Wc1_1 Wc1_2]);

for k = 1:num_seg
    x = icp_seg(k,:);
    xlp = filtfilt(B1,A1,x); % LPF

    % ----------------------- PSD estimate ---------------------------------
    % Welch method (takes segments of signal -- heart rate changes over time). 

    r = 10;
    xlpd = decimate(xlp,r);
    fsd = fix(fs/r);
    
    [pxx,fp] = pwelch(xlpd,13,[],2^12,fsd);

    % ---- find maximumum power freq and freqeuncies on either side ------ 
    % Additional frequencies at 85% power
    imax = find(pxx==max(pxx),1);
    fheart_max = fp(imax); % -------------------------------------------------------------------------- NNEED?
    plim = 0.85*max(pxx);
    
    % split data to find frequencies
    pxx1 = pxx(1:imax);
    pxx2 = pxx(imax+1:end);
    
    i1 = find(pxx1>=plim,1); % find first point at or above 85% (data increasing)
    fh1 = fp(i1);
    i2 = find(pxx2<=plim,1); % find first point at or below 85% (data decreasing)
    i2 = i2 + imax;          % since data was split, add halfway point to be correct
    fh2 = fp(i2);
    
    % ------ Notch filter (keep cardiac component @ max power freq) --------
    
    Wc2_1 = fh1/(fs/2);
    Wc2_2 = fh2/(fs/2);
    [B2,A2] = ellip(2,0.5,40,[Wc2_1, Wc2_2]); % tested stability for range of freqs (ok)
    
    xnotch = filtfilt(B2,A2,x);

    icp_seg(k,:) = xnotch; % replace with processed segment
end 

% --------------- Make matrix icp_seg back into array ------------------

icp_filt_full = reshape(icp_seg.',1,[]);


% -------- ----------------------derivative -----------------------------

icp_dv = [diff(icp_filt_full), 0]; % append zero to deriv to make same length as icp 

% ------------ Moving average filter to smooth segment transitions -----

icp_smooth = icp_dv;
M = 34;
B3 = ones(1,M)/M;

for k = 1:num_seg-1
    seg_end = k*seg_len-fix(0.42*fs);       % ~ 0.42 sec BEFORE segment end 
    next_seg_strt = k*seg_len+fix(0.42*fs); % ~ 0.42 sec AFTER next segment start 
    icp_smooth_seg = filtfilt(B3,1,icp_dv(seg_end:next_seg_strt));  % moving average filter for the segment transitions
    icp_smooth(seg_end:next_seg_strt) = icp_smooth_seg;  % replace abrupt transition with smoothed one
end 

figure('Color',[1 1 1]);
h = plot(t,icp_dv);
set(h,'LineWidth',1.2);
% set(h,'Color',[0 0.6 0.6]);
hold on;
h = plot(t,icp_smooth);
set(h,'LineWidth',1.2);
% set(h,'Color',[0 0.6 0.6]);
xlabel('Time (s)');
ylabel('ICP (mmHG)');
title('Segment Transitions');
legend('Raw derivative', 'Smoothed signal');

% -------- Detect points in between max and min of derivative ----------

a = icp_smooth;
L = length(icp_smooth);
idmax = find((a(1:L-2)<a(2:L-1))&(a(2:L-1)>a(3:L)))+1;
idmin = find((a(1:L-2)>a(2:L-1))&(a(2:L-1)<a(3:L)))+1;

ppi = zeros(1,numel(idmax));

for k = 1:numel(idmax)
    i_start = idmax(k);

    if k < numel(idmin)  % if the interval end is greater than the length of signal, just use the end of signal
        i_end = idmin(k+1);
    else
        i_end = numel(icp_smooth);
    end

    x_intrvl = sig(i_start:i_end);
    pp = find(x_intrvl==max(x_intrvl),1);
    ppi(k) = pp+i_start-1;
end

figure('Color',[1 1 1]);
h = plot(t,sig);
set(h,'LineWidth',1.2);
set(h,'Color',[0.8 0.8 1]);
hold on;

h = plot(t,icp_filt_full+mean(sig)-0.5);
set(h,'LineWidth',1.2);
set(h,'Color',[1 0.6 0.6]);
hold on;

h = plot(t,icp_smooth+mean(sig)-0.5);
set(h,'LineWidth',1.2);
set(h,'Color',[0 0.6 0.6]);
hold on;

h = plot(t,icp_dv+mean(sig)-0.5);
set(h,'LineWidth',1.2);
set(h,'Color',[0.2 0.73 1]);
hold on;

h = plot(idmax/fs,a(idmax)+mean(sig)-0.5,'r.');
set(h,'MarkerSize',18);
hold on;
h = plot(idmin/fs,a(idmin)+mean(sig)-0.5,'b.');
set(h,'MarkerSize',18);
hold on;

h = plot(ppi/fs,sig(ppi),'g.');
set(h,'MarkerSize',18); 
legend('x','xfilt','dvsmooth','dv')


figure('Color',[1 1 1]);
h = plot(t,sig);
set(h,'LineWidth',1.2);
set(h,'Color',[0.8 0.8 1]);
hold on;
h = plot(ppi/fs,sig(ppi),'g.');
set(h,'MarkerSize',18); 


%%
figure('Color',[1 1 1]);
h = plot(t,sig);          % signal
set(h,'LineWidth',1.2);
set(h,'Color',[0.8 0.8 1]);
hold on;

h = plot(t,icp_filt_full+mean(sig)-0.5);  % Cardic Component
set(h,'LineWidth',1.2);
set(h,'Color',[.20 0.47 1]);
hold on;

h = plot(t,icp_smooth+mean(sig)-0.5);  % Smoothed transition
set(h,'LineWidth',1.2);
set(h,'Color',[1 0.2 1]);
hold on;

h = plot(idmax/fs,a(idmax)+mean(sig)-0.5,'g.');
set(h,'MarkerSize',18);
hold on;
h = plot(idmin/fs,a(idmin)+mean(sig)-0.5,'r.');
set(h,'MarkerSize',18);
hold on;

h = plot(ppi/fs,sig(ppi),'b*');
set(h,'MarkerSize',5); 
xlim([9.6 11]);
xlabel('Time (s)');
ylabel('ICP (mmHG)');
legend('ICP','Cardiac Component','Smoothed Transition','Max of DV','Min of DV', 'Percussion Peak')
