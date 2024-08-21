%% Introduction to Programming
% by Dominic Meads
%
%% Common Elements to High-Level Programming Languages
% High-level: with abstraction higher than assembly
%
% * Variables, Data Types, Data Structures
% * Input/Output (Get data from keyboard, file device; output to screen)
% * Mathematical & Logical Operators
% * Conditional Execution (if .. elseif .. else)
% * Repetition (Iteration) / Loops (for or while)
% * Functional Programming
% * Objects
% * Native & Specializeed Libraries (Toolboxes)
% * GUIs

%% Variables
% A variable is a storage location paired with an associated symbolic
% name (identifier) that contains some quantity of information (value)
%
% It is helpful to think of a variable as a container in a program that 
% holds information, that is, a container in which a data value can be
% stored inside the computer's memory. 
%
% Variable Concepts
%
% * Creating Variables
% * Data Types (e.g., single, double, char, logical)
% * Data Structures (e.g., scalars, chars. 1D/2D arrats, vector, matrix)
% * MATLAB Constants (e.g., pi, inf, i, j)
% * Bulit-in Functions (MATLAB) to work with variables (who, whos, clear)
%
clc

% Scalar
a1 = 4

% 1D Array = Vector
a2 = [1 2 3 4 5]

% Character
a3 = 'a'

% Vector Using Colon Operator
a4 = 1:5

% Transpose
a5 = a4'

% String (Vector 1D, 1 Row Matrix)
a6 = 'Dominic'

% Accessing Elements
a7 = a6(2)
a8 = a6(1:3)

% Creating a Matrix
a9 = [1 2 3; 4 5 6; 8 9 10]

% Accessing Elements in a Matrix
a10 = a9(1:2, end)

% Constants
a11 = pi
a12 = j
a13 = NaN
a14 = inf

% Data Types and Conversion of Data Types
a15 = single(a1)

% Special Variables * Functions
a16 = []
a17 = zeros(4,4)
a18 = ones(4,4)
a19 = diag(1:4)

% Replace Variable, Segments of Variables
a20 = a19;
a20(end, :) = ones(1,4)
a20(:, end) = 1:4
a20(2:3, 2:3) = [1, 2; 3, 4]

% Special Functions to Work With Variables
length(a2)
size(a2)
who
whos
%clear

%% Input/Output
% Input from keyboard, file, device &
% Output to screen, file, device
%
% Key Functions and Concepts
%
% * Input
% * Disp
% * fprintf, fopen, fclose
% * format
% * load
%

% Input from Screen
x = 5;
y = 6.6;
name = 'Dominic';


%x = input('Enter the value of x: ')
%y = input('Enter the value of y: ')
%name = input('Enter your name: ', 's')
z = x+y;

% Simple Screen Output
disp('The sum of the variables is: ')
z

disp(['The sum of the variables is: ', num2str(z)])

% Output to Screen (More powerful than disp)
fprintf('The sum of the variables is %g %s \n', z, name)

% Output to File
fid = fopen('myResults.txt', 'wt');
fprintf(fid, 'The sum of the variables is %g %s \n', z, name);
fclose(fid);

% Format
format short
format long
format bank
format short e

% Read data from a file
load handel
save('MyData', 'x', 'y', 'z')

%% Mathematical and logical operators
% High level languages all contain operators:
% assignment, arithmeticr, relational, and logical

% Arithmetic Operators
a1 = 5+5
a2 = 5-5
a3 = 5*5
a4 = 5/5

% Arithmetic Vector (element by element)
a6 = 1:5
a7 = 1:5
a8 = a6+a7
a9 = a6-a7
a10 = a6.*a7
a11 = a6./a7
a12 = a6.^a7

% Arithmetic Matrix
a13 = a6*a7'
a14 = a6'*a7
a15 = ones(3,3)*ones(3,3)

% Relational (Comparison) Operators
a = rand(1,5)
b = rand(1,5)

x_r1 = a < b
x_r2 = a > b
x_r3 = a <= b
x_r4 = a >= b
x_r5 = a == b 
x_r6 = a ~= b

% Logical Operators

a = [0 1 0 1]
b = [1 0 1 1]

z_and = a & b
z_or = a | b
z_not = ~a

%% Conditional Execution
% Relational and Logical Operators are combined with an "if..else.."
% statement structure to selectively execute blocks of code by performing 
% a basic conditonal test that evaluates a given
% expression for a boolean value of true (1) or false (0)

% if statement

if 1 > 5
    disp('Yes, 5 is greater than 1')
    disp('Thanks for asking')
end

% if else statement

a = rand
b = rand

if a < b 
    disp('a is less than b')
    a
    b
else
    disp('a is greater than b')
    a
    b
end 

% if elseif elseif elseif... else

a = rand
b = rand

if a < b 
    disp('a is less than b')
    a
    b
elseif a > b
    disp('a is greater than b')
    a
    b
elseif a == b
    disp('a is equal to than b')
    a
    b 
else
    disp('error')
end 


%% Repetition and Iteration
% Repeated evalutation of a block of code (a set of statements) is
% acomplished using "for loops" and "while loops"
%
% Key concepts
% 
% * For loops (determined # of repetitions)
% * While loops (undetermined # of repetitions) 
% * Vectorization (avoid loops by using faster vector operations)

% For loop

clc
fid = fopen('myResults1.txt', 'wt');
for k = 1:10
    disp('Program entered the for loop');
    fprintf('For loop iteration %d \n', k);
    fprintf(fid,'For loop iteration %d \n', k);
end
fclose(fid);

% While loop
fid = fopen('myResults2.txt', 'wt');
k = 1
while k <= 10
    disp('Program entered the while loop');
    fprintf('while loop iteration %d \n', k);
    fprintf(fid,'while loop iteration %d \n', k);
    k = k+1;
end
fclose(fid);

% While loop with unkown no of iterations

% fid = fopen('myResults3.txt', 'wt');
% k1 = 1;
% k = 1;
% while k1 ~= 0
%     disp('Program entered the while loop');
%     fprintf('while loop iteration %d \n', k);
%     fprintf(fid,'while loop iteration %d \n', k);
%     k1 = input('Do you want to stop [Enter 0]')
%     k = k+1;
% end
% fclose(fid);

% For loops with noninteger increments

y = zeros(1,16); % pre-allocate memory
k = 1;
for x = 0:pi/15:pi
    y(k) = cos(x)
    k = k+1;
end

% Vectorization (no for loops needed)

x = 0:pi/15:pi
y = cos(x)

%% User-defined functions
% Functions are used to enclose a section of code that provides specific
% functionality to the program (sub programs). User-defines functions work
% just like built in functions
%
% Advantages o fusing functions include:
%
% * Addition specific functionality that can be re-used
% * Make programs easier to understand and mantain
% * Ability to divide the work of larger projects among programmers by
%   working on different function with well defined interfaces
%
% Key Concepts
%
% * Declaration of User Define Functions
% * Definition of User define Functions
% * Variable scope
% * Function arguments (input parameters) & passing Arguments by value
% * Output Parameters

a = [1 2 3; 4 5 6; 7 8 9]
b = ones(3,3)

[sum,product] = addmult(a,b)

%% Reading from file

fid = fopen('USGS_discharge_RR_Mcleod.txt','r+');
header = 28;  % 28 header lines
data = textscan(fid,'%s%s%s%s%s%d%s', header, 'Delimiter','    ');
discharge = data{6};
fclose(fid);
figure;
h = plot(discharge);
ylabel('Discharge (cubic-ft/sec)');
