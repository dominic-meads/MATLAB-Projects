function fi = PressureDetect(x,fs,pf)
% PressureDetect: Returns the sample numbers (indecies) of the percussion 
%                 peak of an ICP xnal. 
%
% function fi = PressureDetect(x,fs,pf)
%
% INPUTS:
% x  = ICP signal
% fs = Sampling frequency (default 125 Hz)
% pf = Plot flag: 0 = no plot (default), 1 = plot to screen
%
% OUTPUTS:
% fi = Percussion peak (P1) index, samples
%
% EXAMPLE:
% load ICP.mat;
% x = icp;
% fs = 125;
% fi = PressureDetect(x,fs,1);
% 
% Returns percussion peak indecies and optionally plots signal/peaks to screen
%
% Version 1.0 by Dominic Meads

    % ------------------- handel default inputs and empty arguments -------
    % https://www.mathworks.com/matlabcentral/answers/1951763

    arguments
        x
        fs 
        pf
    end

    if nargin < 2 || isempty(fs)
        fs = 125;
    end

    if nargin < 3 || isempty(pf)
        pf = 0;
    end

    len = length(x);
    t = (0:length(x)-1)/fs;
    win_time = 10;
    seg_len = win_time*fs;
    num_seg = ceil(length(x)/seg_len);
    icp_seg = zeros(num_seg,seg_len);
    
    % --------------------- fill matrix with icp xnal --------------------
    
    for k = 1:num_seg
        i_start = (k-1)*seg_len+1;               % start index
        i_stop = k*seg_len;                      % stop index
        if i_stop <= len                         % fill matrix with the entire segment
            icp_seg(k,:) = x(i_start:i_stop);
        else                                     % fill matrix with remainder of xnal
            rem = len-((k-1)*seg_len);
            icp_seg(k,1:rem) = x(i_start:len);
        end 
    end
    
    % ----------------- Filter low freq trend and high freq noise ---------- 
    % Min heart rate 27 bpm (0.45 Hz), fc1_1 = 0.45 
    % link: World record (google)
    % Max heart rate (Ventricular Tachycardia) 250 bpm (4.17 Hz), fc1_2 = 4.17
    % link: https://ecgwaves.com/topic/ventricular-tachycardia-vt-ecg-treatment-causes-management/
    
    fc1_1 = 0.45;
    fc1_2 = 4.17;
    Wc1_1 = fc1_1/(fs/2);
    Wc1_2 = fc1_2/(fs/2);
    [B1,A1] = ellip(4,0.5,40,[Wc1_1 Wc1_2]);
    
    for k = 1:num_seg
        x1 = icp_seg(k,:);
        xlp = filtfilt(B1,A1,x1); % LPF
    
        % ----------------------- PSD estimate ---------------------------------
        % Welch method (takes segments of xnal -- heart rate changes over time). 
    
        r = 10;
        xlpd = decimate(xlp,r);
        fsd = fix(fs/r);
        
        [pxx,fp] = pwelch(xlpd,13,[],2^12,fsd);
    
        % ---- find maximumum power freq and freqeuncies on either side ------ 
        % Additional frequencies at 85% power
        imax = find(pxx==max(pxx),1);
        plim = 0.85*max(pxx);
        
        % split data to find frequencies
        pxx1 = pxx(1:imax);
        pxx2 = pxx(imax+1:end);
        
        i1 = find(pxx1>=plim,1); % find first point at or above 85% (data increasing)
        fh1 = fp(i1);
        i2 = find(pxx2<=plim,1); % find first point at or below 85% (data decreasing)
        i2 = i2 + imax;          % since data was split, add halfway point to be correct
        fh2 = fp(i2);
        
        % ------ Peak filter (keep cardiac component @ max power freq) --------
        
        Wc2_1 = fh1/(fs/2);
        Wc2_2 = fh2/(fs/2);
        [B2,A2] = ellip(2,0.5,40,[Wc2_1, Wc2_2]); % tested stability for range of freqs (ok)
        
        xpeak = filtfilt(B2,A2,x1);
    
        icp_seg(k,:) = xpeak; % replace with processed segment
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
    
    % -------- Detect points in between max and min of derivative ----------
    
    a = icp_smooth;
    L = length(icp_smooth);
    idmax = find((a(1:L-2)<a(2:L-1))&(a(2:L-1)>a(3:L)))+1;
    idmin = find((a(1:L-2)>a(2:L-1))&(a(2:L-1)<a(3:L)))+1;
    
    fi = zeros(1,numel(idmax));
    
    for k = 1:numel(idmax)
        i_start = idmax(k);
    
        if k < numel(idmin)  % if the interval end is greater than the length of signal, just use the end of signal
            i_end = idmin(k+1);
        else
            i_end = numel(icp_smooth);
        end
    
        x_intrvl = x(i_start:i_end);
        pp = find(x_intrvl==max(x_intrvl),1);
        fi(k) = pp+i_start-1;
    end
    
    % ------------------------- Determine plotting with pf ----------------
    if pf == 1
        figure('Color',[1 1 1]);
        h = plot(t,x);
        set(h,'LineWidth',1.2);
        set(h,'Color',[0.2 0.73 1]);
        hold on;
        h = plot(fi/fs,x(fi),'r.');
        set(h,'MarkerSize',15);
        title('ICP Signal with Location of Percussion Peaks');
        legend('ICP', 'P1 location');
        ylabel('ICP (mmHG)');
        xlabel('Time (s)');
        xlim([min(t) max(t)]);
    end


end 


   
