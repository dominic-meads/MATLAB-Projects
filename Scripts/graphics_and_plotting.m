%% Lab 2 
% by Dominic Meads

close all
clear all
clc

%% 2). Different plots of functions

x = 1:0.1:10;
y1 = x;
y2 = x.^2;
y3 = x.^3;
y4 = sqrt(x);
y5 = exp(x);
y6 = log10(x);

figure('Color',[1 1 1]);
plot(x,y1);
figure('Color',[1 1 1]);
plot(x,y2);
figure('Color',[1 1 1]);
plot(x,y3);
figure('Color',[1 1 1]);
plot(x,y4);
figure('Color',[1 1 1]);
plot(x,y5);
figure('Color',[1 1 1]);
plot(x,y6);

%% 3). plot x(t) = 10cos(2pi5t) + 0.5cos(2pi60t)

fmax = 60;
fs = 12*fmax; 
Ts = 1/fs; 
n = 0:1*fs; % sample for 1 second
t = n*Ts;
x = 10*cos(2*pi*5*t)+0.5*cos(2*pi*60*t);

figure('Color',[1 1 1]);
h = plot(t,x);
set(h,'Linewidth', 1);
set(h,'Color',[0.4 0.4 1]);
hold on;
h = plot(t,x,'.');
xlabel('Time (s)');
ylabel('Amplitude');
set(h,'MarkerSize', 10);
set(h,'Color',[1 0.3 0.3]);
print -dtiff -r300 Lab2_sinusoid; 

fmax = 60;
fs = 2*fmax; % reduce sampling frequency
Ts = 1/fs; 
n = 0:1*fs; 
t = n*Ts;
x = 10*cos(2*pi*5*t)+0.5*cos(2*pi*60*t);

figure('Color',[1 1 1]);
h = plot(t,x);
set(h,'Linewidth', 1);
set(h,'Color',[0.4 0.4 1]);
hold on;
h = plot(t,x,'.');
title('Graph of x(t) with reduced sampling frequency')
xlabel('Time (s)');
ylabel('Amplitude');
set(h,'MarkerSize', 10);
set(h,'Color',[1 0.3 0.3]);

%% 4). Create 6 subplots of x(t)

figure('Color',[1 1 1]);

subplot(3,2,1);
h = plot(x);
set(h,'Color',[0.4 1 0.4]);
grid off;

subplot(3,2,2);
h = plot(x);
set(h,'Linewidth', 4);
grid on;

subplot(3,2,3);
h = plot(x);
set(h,'LineStyle', '--');
axis tight;

subplot(3,2,4);
stem(x);
axis tight;

subplot(3,2,5);
h = stem(x);
set(h,'MarkerSize', 4);
set(h,'Color',[1 0.3 0.3]);

subplot(3,2,6);
plot(x);
hold on;
h = plot(x,'.');
set(h,'Marker', '+');

%% 5). Plot the function f(x)

x1 = 0:0.2:25;
f = (1./(((x1-3).^2)+0.01))+(1./(((x1-0.9).^2)+0.04))-6;

figure('Color',[1 1 1]);
plot(x1,f);

%% 6). Subplots

V = 0:0.1:10;
R = 1000;
I = V/R;

figure('Color',[1 1 1]);
subplot(2,2,1);
h1 = plot(V,I);
title('I-V Characteristics of 1KOhm Resistor');
xlabel('Voltage (V)');
ylabel('Current (A)');

subplot(2,2,2)
ysquare = x.^2;
h2 = plot(t,ysquare);
title('x^2(t)');

subplot(2,2,3);
yabs = abs(x);
h3 = plot(t,yabs);
title('|x(t)|');

subplot(2,2,4);
yexp = exp(-5*t).*x;
h4 = plot(t,yexp);
title('e^{-5t} * x(t)');

%% 7). Use get, set, and gca

hold on;
get(h4)
set(h4,'LineStyle','--');
gca % properties of axes

%% 8). Normal random numbers

r1 = randn(1,1000);

figure('Color',[1 1 1]);
subplot(1,2,1);
plot(r1);
title('Time Series');
subplot(1,2,2);
hist(r1);
title('Histogram');

%% 9). Uniform random numbers

r2 = rand(1,1000);

figure('Color',[1 1 1]);
subplot(1,2,1);
plot(r2);
title('Time Series');
subplot(1,2,2);
hist(r2);
title('Histogram');

%% 10). Uniform random numbers part 2

r3 = rand(10,1000);

figure('Color',[1 1 1]);
subplot(5,2,1);
hist(r3(1,:));
subplot(5,2,2);
hist(r3(2,:));
subplot(5,2,3);
hist(r3(3,:));
subplot(5,2,4);
hist(r3(4,:));
subplot(5,2,5);
hist(r3(5,:));
subplot(5,2,6);
hist(r3(6,:));
subplot(5,2,7);
hist(r3(7,:));
subplot(5,2,8);
hist(r3(8,:));
subplot(5,2,9);
hist(r3(9,:));
subplot(5,2,10);
hist(r3(10,:));

figure('Color', [1 1 1]);
hist(sum(r3));
title('Sum of 10 Instances of Uniformly Distributed Random Numbers');

%% 11). Load Data

fid = fopen('USGS_discharge_RR_Mcleod.txt','r');

%% 3D Graphics

% 1). Plotting

[X,Y] = meshgrid(-2.5:.2:2.5, -2.5:.2:2.5);
Z = (X.^3 + 3.*Y^3).*exp(1 - X.^2 - Y.^2);
figure('Color', [1 1 1]);
mesh(X,Y,Z);

% 2). Viewpoint changes
figure('Color', [1 1 1]);
mesh(X,Y,Z);
view(-15,60);
rotate3d on;

% 3). contour plot
figure('Color', [1 1 1]);
contour(X,Y,Z);

% 4). Different types of plots
figure('Color', [1 1 1]);
subplot(2,2,1);
surf(X,Y,Z);
title('surf');
subplot(2,2,2);
meshc(X,Y,Z);
title('meshc');
subplot(2,2,3);
meshz(X,Y,Z);
title('meshz');
subplot(2,2,4);
waterfall(X,Y,Z);
title('waterfall');

% 5). surfl
figure('Color', [1 1 1]);
surfl(X,Y,Z);
colorbar;

% 6). MATLAB Logo
logo;
cameramenu


%% Solving Math Problems

% 1). Euler's Equation
ans = exp(j*pi)+1
abs(ans) % very close to zero

% 2). Euler's Formula

theta = 0:pi/8:pi;
arr1 = exp(j.*theta);
arr2 = cos(theta) + j*sin(theta);
if isequal(arr1,arr2) == 1
    disp('Eulers formula is proved');
else
    disp('Error');
end

% EXTRA CREDIT showing plots are equal 
figure('Color',[1 1 1]);
angle = 0:0.1:4*2*pi;  % four revolutions
x1 = exp(i*angle);
x2 = cos(angle)+j*sin(angle);
subplot(2,2,1);
h = plot(real(x1));
title('Re[e^j*theta]');
ylabel('Real(x)');
axis tight;
subplot(2,2,2);
h = plot(imag(x1));
title('Im[e^j*theta]');
ylabel('imag(x)');
axis tight;
subplot(2,2,3);
h = plot(real(x2));
title('Re[cos(theta)+jsin(theta)]');
ylabel('real(x)');
axis tight;
subplot(2,2,4);
h = plot(imag(x2));
title('Im[cos(theta)+jsin(theta)]');
ylabel('Imag(x)');
axis tight;


% 3). Integrals
dx = 0.0001;
x = 0:dx:1;

I1 = sum(x*dx)
disp('answer should be ~0.5')
I2 = sum(x.^2*dx)
disp('answer should be ~0.3333')
I3 = sum(x.^3*dx)
disp('answer should be ~0.25')
I4 = sum(exp(x)*dx)
disp('answer should be ~1.71828')


% 4). More integrals
y = 1./(sqrt(x.^2+2));
I5 = sum(y*dx)
disp('answer should be ~0.65848')
figure;
plot(x,y);
axis tight;

% 5). 

dx = 0.01;
x = 0:dx:2*pi;
y = abs(sin(x));
I = find(y > 0.9);
length(I)*dx/(2*pi)
