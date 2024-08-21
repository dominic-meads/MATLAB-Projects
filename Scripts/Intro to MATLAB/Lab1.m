%% Lab 1 Basic MATLAB functions and Matrix Maniuplations
% Dominic Meads

%% 1). & 2). Script with MATLAB Functions
clc
clear all
a = 1;
b = 2;
c = 3;
who
demo % bring up demo
D = ones(2,2) % 2x2 matrix with ones only
E = zeros(2,2) % 2x2 matrix with zeros only
F = eye(2,2) % 2x2 matrix with ones on diagonal
H = diag(1:5) 

Result = input('What is 2+2? ');
pause(3);
disp('If you entered 4, you are correct!')

fid = fopen('myFile.txt','wt');
fprintf(fid,'my name is Dom');
fclose(fid)
save Lab1
clear all
clc

load Lab1

format long
pi
x = 5.364554773737333;
format short
x

x = 0:0.1:10;
y = x.^3;
plot(y);
xlabel('x axis');
ylabel('y axis');
title('my graph');

length(x)
size(x)

a = [1 2 3 1 1 3 2 1 4];
me = mean(a)
md = median(a)
stdev = std(a)
ma = max(a)
mi = min(a)
s = sum(a)

n = -4
abs(n)
b = sqrt(abs(n))
d = exp(2)
cosine = cos(0)
sine = sin(0)
poly = [1 0 -9];
r = roots(poly)

sort(a)

figure('Color',[1 1 1]);
r = rand(1,1000);
hist(r)

any(r)
r2 = zeros(10,1)
any(r2)

p = [0 0 0 1 0 2]
I = find(p)

m = ones(1,5)
all(m)

v = [1 2 3 4 5 6]
vdiff = diff(v)

%% 3). Creating Matricies

% Create matrix A using only 3 lines of code.
% Use only diag function and colon/transpose operators 
A = diag(1:10);  % create 10x10 matrix with 1:10 on diagonal
A(1,:) = [1:10]; % replace first row (all columns) with 1:10
A(:,1) = A(1,:)' % replace all columns in column 1 w/transposed first row 

% create matrix M using only 1 line of code
% Use only eye, ones, and zeros functions
M = [2*eye(2), zeros(2,3);3*ones(2,3),4*ones(2,2)]

% Explanation: create two matricies and combine using semi-colon
M1 = [2*eye(2), zeros(2,3)]
M2 = [3*ones(2,3),4*ones(2,2)] 
m = [M1;M2]

%% 4). Computing Factorials
% compute using function cumprod and single line of code

factorial_list = [1:10;cumprod(1:10)]'

% explanation: similar to problem 3, create two vectors and combine.
% transpose matrix to display as shown in the problem.

%% Conceptual Questions

% a). MATLAB costs money (octave is free), however, MATLAB has some
%     built-in functions and toolboxes that Octave does not have.

% -------------------------------------------------------------------

% b). MATLAB functions are built into the development environment. A
%     script is user-created and uses these functions to do something.

% -------------------------------------------------------------------

% c). clc just clears the command window. If I have variables
a = 1
b = 2
c = 3
% then run it, the result will be displayed in the command window. 
% clc will clear this command window, however, the variables will remain in
% memory. 
clc
% clear will clear the variables from memory.
clear 
% the variables no longer exist

% ----------------------------------------------------------------------

% d). I will create two arrays:
arr1 = [1 2 3 4];
arr2 = [2 2 2 2];
% '*' will treat the arrays as matricies, and when I try to multiply, it
% will not work (I must transpose arr2 to satisfy rules of matrix
% multiplication). This calculates a dot product
arr1*arr2'
% additionally, '*' can be used for multiplying single variables (1x1
% matrix)
k = 4;
l = 5;
k*l
% '.*' is element by element multiplication, meanining it will multiply
% each induvidial element by each other. arr1(1)*arr2(1), arr1(2)*arr2(2),
% etc.
arr1.*arr2




