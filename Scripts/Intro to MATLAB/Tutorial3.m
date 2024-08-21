%% Tutorial 3
% by Dominic Meads

%% Using the sum function to approximate integrals
% definite integrals

dx = 0.0001; % differnetial intervals
x = 0:dx:1;
y = x;
A1 = sum(y*dx) % close to actual answer
y = x.^2;
A2 = sum(y*dx)

%% Plot of Euler Equation
% Euler

figure('Color',[1 1 1]);
angle = 0:0.1:4*2*pi;  % four revolutions
x = exp(i*angle);
subplot(3,1,1);
h = plot(x);
xlabel('Real(x)');
ylabel('Imag(x)');
subplot(3,1,2);
h = plot(real(x));
ylabel('Real(x)');
axis tight;
subplot(3,1,3);
h = plot(imag(x));
ylabel('Imag(x)');
axis tight;

figure('Color',[1 1 1]);
h = plot3(angle,real(x),imag(x));
xlabel('Angle');
ylabel('Real(x)');
zlabel('Imag(x)');


%% Introduction to For Loops
% Visualization of the Euler Equation

figure('Color',[1 1 1]);
for k = 1:length(angle);
    h = plot3(angle(k),real(x(k)),imag(x(k)),'.');
    xlabel('Angle');
    ylabel('Real(x)');
    zlabel('Imag(x)');
    hold on;
    pause(0.001);
end;













