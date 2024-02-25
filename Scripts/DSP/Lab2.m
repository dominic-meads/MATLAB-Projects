%% Lab 2 
% by Dominic Meads
clc
close all
clear 
%% x(n) = a^n

n = -10:10;
a = [2 -2 -0.5 0.5];
x = zeros(4,length(n));
figure('Color', [1 1 1]);
for k = 1:length(a)
    x(k,:) = a(k).^n;
    subplot(2,2,k);
    plot(n,x(k,:));
    title(['x(n) = ',num2str(a(k)),'^n']);
    xlabel('n');
    ylabel('x(n)')
end 

% When a is postive, we see a similar shape to an
% exponential increase (a > 1), and exponential decay (a < 1). When a is
% negative, odd n values result in a negative x(n), while even n values
% result in a postive x(n). 

%% x(n) = e^(alpha + j*beta)n

alpha = [1.1 -1.1 0.9 -0.9];
beta = [1.1 -1.1 0.9 -0.9];
x = zeros(4,length(n));

% Real parts
figure('Color', [1 1 1]);
for k = 1:length(alpha)
    x(k,:) = exp((alpha(k)+j*beta(k)).*n);
    subplot(2,2,k);
    plot(n,real(x(k,:)));
    title(['Re[e^{(',num2str(alpha(k)),'+j*',num2str(beta(k)),')n}]']);
    xlabel('n');
    ylabel('x(n)');
end 

% Imaginary parts
figure('Color', [1 1 1]);
for k = 1:length(alpha)
    x(k,:) = exp((alpha(k)+j*beta(k)).*n);
    subplot(2,2,k);
    plot(n,imag(x(k,:)));
    title(['Im[e^{(',num2str(alpha(k)),'+j*',num2str(beta(k)),')n}]']);
    xlabel('n');
    ylabel('x(n)');
end 

% Magnitude
figure('Color', [1 1 1]);
for k = 1:length(alpha)
    x(k,:) = exp((alpha(k)+j*beta(k)).*n);
    subplot(2,2,k);
    plot(n,abs(x(k,:)));
    title(['Magnitude[e^{(',num2str(alpha(k)),'+j*',num2str(beta(k)),')n}]']);
    xlabel('n');
    ylabel('|x(n)|');
end 

% Angle
figure('Color', [1 1 1]);
for k = 1:length(alpha)
    x(k,:) = exp((alpha(k)+j*beta(k)).*n);
    subplot(2,2,k);
    plot(n,angle(x(k,:)));
    title(['Angle[e^{(',num2str(alpha(k)),'+j*',num2str(beta(k)),')n}]']);
    xlabel('n');
    ylabel('Angle[x(n)] (rad)');
end 

% I chose to explore four cases where alpha and beta were both postive and
% greater than one, both negative and greater than 1, both positive and
% less than 1, and both negative and less than 1. The postitive and
% negative alphas/beta graphs for both greater than 1 and less than 1 were
% mirrored. The real parts of each are negative for the extreme parts of n,
% while in the imaginary parts, x(n) is negative for extreme n when
% alpha/beta > 1. x(n) is positive for extreme n when alpha/beta < 1. All
% magnitude plots are similarly shaped. If we imagine the phase graphs as
% sawtooth waveforms, apha/beta > 1 has a higher "frequency"

%% 