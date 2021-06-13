% [Sascha Spors](https://orcid.org/0000-0001-7225-9992),
% Professorship Signal Theory and Digital Signal Processing,
% [Institute of Communications Engineering (INT)](https://www.int.uni-rostock.de/),
% Faculty of Computer Science and Electrical Engineering (IEF),
% [University of Rostock, Germany](https://www.uni-rostock.de/en/)
% 
% # Tutorial Signals and Systems (Signal- und Systemtheorie)
% 
% Summer Semester 2021 (Bachelor Course #24015)
% 
% - lecture: https://github.com/spatialaudio/signals-and-systems-lecture
% - tutorial: https://github.com/spatialaudio/signals-and-systems-exercises
% 
% WIP...
% The project is currently under heavy development while adding new material for the summer semester 2021
% 
% Feel free to contact lecturer [frank.schultz@uni-rostock.de](https://orcid.org/0000-0002-3010-0294)

% Exercise 8.2

clear all
close all
clc

% signal
x = [1, 1, 2, -1]  % non-zero elements
Nx = size(x,2)
kx = 1  % start index for first non-zero entry

% signal, we can interpret this as impulse response of LTI system
h = [2, 1, -1]
Nh = size(h,2)
kh = 0

% convolution
Ny = Nx+Nh-1
ky = kx+kh
y= conv(x,h)


figure
subplot(1,3,1)
k = kx:kx+Nx-1
stem(k,x)
xlabel('k')
ylabel('x[k]')

subplot(1,3,2)
k = kh:kh+Nh-1
stem(k,h)
xlabel('k')
ylabel('h[k]')

subplot(1,3,3)
k = ky:ky+Ny-1
stem(k,y)
xlabel('k')
ylabel('y[k]=x[k]*h[k]')

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
% ``main file, github URL, commit number and/or version tag, year``.
