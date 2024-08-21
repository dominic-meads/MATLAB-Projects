%% Tutorial 4
% by Dominic Meads

close all
clear all
clc
%% Chaos Example
% Chaotic simulator
% Use the for loop to simulate a chaotic system
% defined by the following recursion
%  x(i) = l*x(i-1)(1-x(i-1))
% 3.8 < l < 4.0,
% 0 < x(0) < 1
%
% x = chaos(x0, lambda, vectorLength)

x0 = 0.5;
lambda = 3.8;
vectorLength = 1500;

x = zeros(vectorLength,1);
x(1) = x0;

for k = 2:vectorLength
    x(k) = lambda*x(k-1)*(1-x(k-1));
end

%% Visualizing Chaotic Time-Series
% 2D plots

figure('Color',[1 1 1]);
h = plot(x);
xlabel('Sample Number');
ylabel('Chaotic Number');
box off;

figure('Color',[1 1 1]);
hist(x);

%% Visualizing Chaotic Time-Series part 2
% 2D plots

T = 40;
x1 = x(1:end-T);
x2 = x(T+1:end);

figure('Color',[1 1 1]);
h = plot(x1,x2);
xlabel('x(t)');
ylabel('x(t+T)');

%% Visualizing Chaotic Time-Series part 3
% 3D plots

T = 2;
x1 = x(1:end-2*T);
x2 = x(T+1:end-T);
x3 = x(2*T+1:end);

figure('Color',[1 1 1]);
h = plot3(x1,x2,x3);
xlabel('x(t)');
ylabel('x(t+T)');
zlabel('x(t+2T)');


%% Visualizing Chaotic Time-Series part 4
% 3D plots 

T = 2;
x1 = x(1:end-2*T);
x2 = x(T+1:end-T);
x3 = x(2*T+1:end);

figure('Color',[1 1 1]);

for k = 1:length(x1)
    h = plot3(x1(k),x2(k),x3(k),'.'); hold on;
    set(h,'Markersize', 18);
    xlabel('x(t)');
    ylabel('x(t+T)');
    zlabel('x(t+2T)');
    pause(0.001);
end


