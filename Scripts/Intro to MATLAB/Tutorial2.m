%% Tutorial 2 - Introduction to MATLAB graphics
% by Dominic Meads

close all; clc; clear all;
%% Basic Plotting
% Functions: plot, box off, axis tight; xlabel, ylabel
% get, set, hold on;

x = 0:0.5:10;
y1 = x;
y2 = x.^2;
y3 = sqrt(x);
y4 = log10(x(2:end)); % cannot start log at zero

figure('Color',[1 1 1]);
h = plot(x,y1); % create plot handel
box off; grid on;
hold on;
set(h,'Linewidth', 3);
set(h,'Color',[0.4 0.4 1]); % set color to light blue
h = plot(x,y1,'.'); % plot dots over interpolations
xlabel('x-axis');
ylabel('y-axis');
set(h,'MarkerSize', 18);
set(h,'Color',[1 0.3 0.3]);
axis tight;

figure('Color', [1 1 1]);
h = plot(x,y2); grid on;
box off;  % plots y as a function of x
hold on;
set(h,'Linewidth', 3);
set(h,'Color',[0.4 0.4 1]);
h = plot(x,y2,'.');
xlabel('x-axis');
ylabel('y-axis');
set(h,'MarkerSize', 18);
set(h,'Color',[1 0.3 0.3]);
axis tight;

figure('Color',[1 1 1]);
h = plot(x,y3); grid on;
set(h,'Linewidth', 3);
set(h,'Color',[0.4 0.4 1]);
box off; 
hold on;
h = plot(x,y3,'.');
xlabel('x-axis');
ylabel('y-axis');
set(h,'MarkerSize', 18);
set(h,'Color',[1 0.3 0.3]);
axis tight;

figure('Color',[1 1 1]);
h = plot(x(2:end),y4); grid on;
set(h,'Linewidth', 3);
set(h,'Color',[0.4 0.4 1]);
box off; 
hold on;
h = plot(x(2:end),y4,'.');
xlabel('x-axis');
ylabel('y-axis');
set(h,'MarkerSize', 18);
set(h,'Color',[1 0.3 0.3]);
axis tight;

 %% Subplots
 % Functions: subplot command

figure('Color',  [1 1 1]);
subplot(2,2,1);
h = plot(x,y1); % create plot handel
box off; grid on;
hold on;
set(h,'Linewidth', 3);
set(h,'Color',[0.4 0.4 1]); % set color to light blue
h = plot(x,y1,'.'); % plot dots over interpolations
xlabel('x-axis');
ylabel('y-axis');
set(h,'MarkerSize', 18);
set(h,'Color',[1 0.3 0.3]);
axis tight;

subplot(2,2,2);
h = plot(x,y2); grid on;
box off;  % plots y as a function of x
hold on;
set(h,'Linewidth', 3);
set(h,'Color',[0.4 0.4 1]);
h = plot(x,y2,'.');
xlabel('x-axis');
ylabel('y-axis');
set(h,'MarkerSize', 18);
set(h,'Color',[1 0.3 0.3]);
axis tight;

subplot(2,2,3);
h = plot(x,y3); grid on;
set(h,'Linewidth', 3);
set(h,'Color',[0.4 0.4 1]);
box off; 
hold on;
h = plot(x,y3,'.');
xlabel('x-axis');
ylabel('y-axis');
set(h,'MarkerSize', 18);
set(h,'Color',[1 0.3 0.3]);
axis tight;

subplot(2,2,4);
h = plot(x(2:end),y4); grid on;
set(h,'Linewidth', 3);
set(h,'Color',[0.4 0.4 1]);
box off; 
hold on;
h = plot(x(2:end),y4,'.');
xlabel('x-axis');
ylabel('y-axis');
set(h,'MarkerSize', 18);
set(h,'Color',[1 0.3 0.3]);
axis tight;

% now plot as 1 column
figure('Color',  [1 1 1]);
subplot(4,1,1);
h = plot(x,y1); % create plot handel
box off; grid on;
hold on;
set(h,'Linewidth', 3);
set(h,'Color',[0.4 0.4 1]); % set color to light blue
h = plot(x,y1,'.'); % plot dots over interpolations
xlabel('x-axis');
ylabel('y-axis');
set(h,'MarkerSize', 18);
set(h,'Color',[1 0.3 0.3]);
axis tight;

subplot(4,1,2);
h = plot(x,y2); grid on;
box off;  % plots y as a function of x
hold on;
set(h,'Linewidth', 3);
set(h,'Color',[0.4 0.4 1]);
h = plot(x,y2,'.');
xlabel('x-axis');
ylabel('y-axis');
set(h,'MarkerSize', 18);
set(h,'Color',[1 0.3 0.3]);
axis tight;

subplot(4,1,3);
h = plot(x,y3); grid on;
set(h,'Linewidth', 3);
set(h,'Color',[0.4 0.4 1]);
box off; 
hold on;
h = plot(x,y3,'.');
xlabel('x-axis');
ylabel('y-axis');
set(h,'MarkerSize', 18);
set(h,'Color',[1 0.3 0.3]);
axis tight;

subplot(4,1,4);
h = plot(x(2:end),y4); grid on;
set(h,'Linewidth', 3);
set(h,'Color',[0.4 0.4 1]);
box off; 
hold on;
h = plot(x(2:end),y4,'.');
xlabel('x-axis');
ylabel('y-axis');
set(h,'MarkerSize', 18);
set(h,'Color',[1 0.3 0.3]);
axis tight;

%% Plotting Sinusoidal Signals
% Topic: Sampling Frequency, period, etc.
% x1(t) = 5cos(2pi3t)
% x2(t) = 5cos(2pi3t) + 1cos(2*pi*60*t)

fmax = 3;
fs = 18*fmax; % sampling frequency (3 is fmax which is 3 hz)
Ts = 1/fs; % sampling period
n = 0:2*fs; % samples (for 2 seconds)
t = n*Ts;  % sampling instances
x1 = 5*cos(2*pi*3*t); % sampled signal

figure('Color', [1 1 1]);
h = plot(t,x1);
box off;
set(h,'Linewidth',3);
set(h,'Color',[0.6 0.6 1]);
hold on;
h = plot(t,x1,'.');
box off;
set(h,'MarkerSize',18);
set(h,'Color',[0.1 0.1 1]);
xlabel('Time (s)');
ylabel('Voltage (mV)');

print -dtiff -r300 mysinusoid; % high resultion figure

fmax = 60;
fs = 18*fmax; % sampling frequency (3 is fmax which is 3 hz)
Ts = 1/fs; % sampling period
n = 0:2*fs; % samples (for 2 seconds)
t = n*Ts;  % sampling instances
x2 = 5*cos(2*pi*3*t) + 1*cos(2*pi*60*t); % sampled signal

figure('Color', [1 1 1]);
h = plot(t,x2);
box off;
set(h,'Linewidth',3);
set(h,'Color',[0.6 0.6 1]);

% dont plot all the markers (makes easier to see)

% hold on;
% h = plot(t,x2,'.');
% box off;
% set(h,'MarkerSize',18);
% set(h,'Color',[0.1 0.1 1]);
% xlabel('Time (s)');
% ylabel('Voltage (mV)');

%% Plotting random numbers, histograms, etc
% Functions: rand, randn; hist

figure('Color', [1 1 1]);
x1 = rand(3,1000);
subplot(3,1,1);
plot(x1(1,:));
subplot(3,1,2);
plot(x1(2,:));
subplot(3,1,3);
plot(x1(3,:));

figure('Color', [1 1 1]);
hist(x1(1,:)); % uniformly distributed

figure('Color', [1 1 1]);
x2 = randn(3,1000); % Gaussian Distribution
subplot(3,1,1);
plot(x2(1,:));
subplot(3,1,2);
plot(x2(2,:));
subplot(3,1,3);
plot(x2(3,:));

figure('Color', [1 1 1]);
hist(x2(1,:)); % higher probability the number is closer to zero

% engineering proof showing that unform noise sources added together are
% normally distributed
figure('Color', [1 1 1]);
hist(sum(x1));

%% 3D plots
% Functions: view, rotate3d, mesh

% 3D MATLAB example
[X,Y] = meshgrid(-2:.2:2, -2:.2:2);
Z = X.*exp(-X.^2 - Y.^2);
figure('Color', [1 1 1]);
surf(X,Y,Z)

figure('Color', [1 1 1]);
mesh(X,Y,Z) % doesnt fill in colors

figure('Color', [1 1 1]);
meshc(X,Y,Z) % projects on to bottom axis

%% MATLAB to solve math problems

exp(i*pi) + 1





































