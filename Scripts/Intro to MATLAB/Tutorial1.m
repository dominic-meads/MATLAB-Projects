clc
A = [1 2 3; 4 5 6;7 8 9]
who
a = A(2,:) % 2nd row, all the columns
b = A(2,2:3) % 2nd row, last two columns
C = A+1;
who
whos
save Tutorial1
clear
load Tutorial1
t = 1:10
t2 = 1:0.5:10 % increment by 0.5
length(t2) % # elements of a vector
size(t2) % size of matrix
cos([0.1 0.01 0.5 0]) % eval cosine at each number in vector
x = 0:0.1:8*pi;
length(x)
y = cos(x);
plot(x,y)
figure
plot(x)
a = [1 2 3]
b = [4 5 6]
c = a.*b %element by element
c = a*b' % dot product
c = sum(a.*b)
a'*b % cross product
figure
plot(x)
y = x.^2;
figure
plot(y)
close all % closes all plots

%% Tutorial 1
% by Dominic Meads
close all

%% Basic MATLAB matrix operations
% This section shows basic matlab operations

A = [1 2 3; 4 5 6;7 8 9]

a = A(2,:)
b = A(:,3)
c = A(1,1:3)

%% Basic plotting
% intro to plot command

t = 0:0.1:4*pi;
y = cos(t);
plot(y)

%% Basic Functions

size(A)
size(b)
length(c)
who
whos

%% Basic Statistics Functions
% Introduction to statistics functions

ma = max(A)
mi = min(A)
u = mean(y)
m = median(y)

x = randn(1000,1);
figure('Color',[1 1 1])
plot(x)
figure('Color',[1 1 1])
hist(x)

mean(x)
std(x)