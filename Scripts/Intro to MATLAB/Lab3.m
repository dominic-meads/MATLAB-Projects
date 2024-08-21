%% Lab 3
% by Dominic Meads

clear all
clc
close all

%% 1).  Function simulation and test

x0 = 0.5;
lambda = 3.8;
vectorLength = 1500;

x = zeros(vectorLength,1);
x(1) = x0;

for k = 2:vectorLength
    x(k) = lambda*x(k-1)*(1-x(k-1));
end

T = 2;
x1 = x(1:end-2*T);
x2 = x(T+1:end-T);
x3 = x(2*T+1:end);

figure('Color',[1 1 1]);
h = plot3(x1,x2,x3);
xlabel('x(t)');
ylabel('x(t+T)');
zlabel('x(t+2T)');

% test
x = chaos(0.9, 3.8, 1000);

%% 2). writing to files

lambda_min = 3.8001; % greater than 3.8
lambda_max = 3.9999; % less than 4

% create array of 30 random lambdas between 3.8 and 4
rand_lambda = (lambda_max-lambda_min).*rand(30,1) + lambda_min;


x0_min = 0.0001;     % greater than zero
x0_max = 0.9999;      % less than 1

% create array of 30 random initial x between 0 and 1
rand_x0 =(x0_max-x0_min).*rand(30,1) + x0_min;

% create directory
mkdir Chaos
Dir = 'C:\Users\demea\OneDrive\Documents\ENGR267\Chaos';

N = 2;
char_N = num2str(N, '%d');
filename = strcat('Chaos',char_N,'.txt');
fid = fopen(fullfile(Dir,filename),'wt');
fprintf(fid,'YOur file has been opened');
fclose(fid);


for k = 1:30
    % run function
    x = chaos(rand_x0(k),rand_lambda(k),1000); 

    % create file name
    char_N = num2str(k, '%d');
    filename = strcat('Chaos',char_N,'.txt');

    % write to file
    Dir = 'C:\Users\demea\OneDrive\Documents\ENGR267\Chaos';
    fid = fopen(fullfile(Dir,filename),'wt');
    fprintf(fid,'%f\n',x);
    fclose(fid);
end

close all;

%% 3). Loading from files and making graphs

x_load = zeros(1000,1);
cd Chaos

for k = 1:30
    % create file name
    char_N = num2str(k, '%d');
    filename = strcat('Chaos',char_N,'.txt');

    % load from file
    x = load(filename);

    % plot and create tiff files
    figure('Color', [1 1 1]);
    plot(x);
    eval(sprintf('print -dtiff chaos%d',k));
    pause(0.1);
end

close all

%% Conceptual Questions

% 1). for loops are best when you have a KNOWN number of iterations (when
% you need to run a program a certain number of times). While loops are
% great when you run a program for an indefinite amount of time (like doing
% something while waiting for input from a keyboard or sensor). 

% 2). Vectorization can be used like a for loop, whitout the for loop
% syntax. it is much more efficient, and uses a vector of values to be
% operated upon. 

% 3). Functionsa are created for usability (something that will be done
% often or something specific). They are called in a script which uses
% different functions to do something.

% 4). Variable scope is kind of like where the variables are in the
% program. You can have "local" variables that only one function may see,
% and you can have global variables that all functions within a script can
% operate on. 


