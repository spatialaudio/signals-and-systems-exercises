% [Sascha Spors](https://orcid.org/0000-0001-7225-9992),
% Professorship Signal Theory and Digital Signal Processing,
% [Institute of Communications Engineering (INT)](https://www.int.uni-rostock.de/),
% Faculty of Computer Science and Electrical Engineering (IEF),
% [University of Rostock, Germany](https://www.uni-rostock.de/en/)
% 
% # Tutorial Signals and Systems (Signal- und Systemtheorie)
% 
% Summer Semester 2022 (Bachelor Course #24015)
% 
% - lecture: https://github.com/spatialaudio/signals-and-systems-lecture
% - tutorial: https://github.com/spatialaudio/signals-and-systems-exercises
% 
% 
% Feel free to contact lecturer [frank.schultz@uni-rostock.de](https://orcid.org/0000-0002-3010-0294)

% Exercise 8.4

clear all
close all
clc

M = 7;
k1 = 0;
x1 = [zeros(1,k1) ones(1,M)];
t1 = [0:length(x1)-1];

N = 4;
k2 = 0;
x2 = [zeros(1,k2) ones(1,N)];
t2 = [0:length(x2)-1];

y = conv(x1,x2,'full');
t = [0:length(y)-1];

L = M+N-1

subplot(3,1,1)
stem(t1,x1, 'r', 'LineWidth',3)
subplot(3,1,2)
stem(t2,x2, 'b')
subplot(3,1,3)
stem(t,y, 'g', 'LineWidth',3)

clc
disp('k1+k2')
k1+k2
disp('k1+k2+y-1')
k1+k2+max(y)-1
disp('k1+k2+L-y')
k1+k2+L-max(y)
disp('k1+k2+L-1')
k1+k2+L-1
disp('|M-N|+1')
abs(M-N)+1
max(y)-1

% ## Copyright
% 
% This tutorial is provided as Open Educational Resource (OER), to be found at
% https://github.com/spatialaudio/signals-and-systems-exercises
% accompanying the OER lecture
% https://github.com/spatialaudio/signals-and-systems-lecture.
% Both are licensed under a) the Creative Commons Attribution 4.0 International
% License for text and graphics and b) the MIT License for source code.
% Please attribute material from the tutorial as *Frank Schultz,
% Continuous- and Discrete-Time Signals and Systems - A Tutorial Featuring
% Computational Examples, University of Rostock* with
% ``github URL, commit number and/or version tag, year, (file name and/or content)``.
