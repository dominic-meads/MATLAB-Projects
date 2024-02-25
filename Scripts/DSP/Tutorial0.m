%% Tutorial 0 - sampling
% by Dominic Meads

%% Plotting a Sinusoid
% x(t) = 5*cos(2*pi*2*t) + cos(2*pi*5*t)

fmax = 2;       % max signal frequency
fs = 3*fmax;   % sampling frequency 
Ts = 1/fs;      % sampling period
n = 0:1*fs-1;     % sampling index
t = n*Ts;       % sampling instances

x = 1*cos(2*pi*2*t);

figure('Color',[1 1 1]);
h = plot(t,x);
hold on;
h = plot(t,x,'r.');
xlabel('Time (s)');
ylabel('Signal');

%% Load ICP and Plot

load ICPcomposite;

figure('Color',[1 1 1]);

t = (0:length(x)-1)/(fs*60);
h = plot(t,x);
xlabel('Time (min)');