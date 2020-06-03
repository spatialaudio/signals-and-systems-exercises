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