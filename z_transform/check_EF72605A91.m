clear all
close all
clc

k = [0:8];

r = 3/4;
phi = pi/3;

a = r*exp(+1j*phi);
xr = a.^k;

a = r*exp(-1j*phi);
xi = a.^k;

x = (xr-xi)/(1j)


subplot(3,1,1)
stem(k, real(xr), 'b'), hold on
stem(k, imag(xr), 'r')
hold off
axis([0,8,-1,1])

subplot(3,1,2)
stem(k, real(xi), 'b'), hold on
stem(k, imag(xi), 'r')
hold off
axis([0,8,-1,1])

subplot(3,1,3)
stem(k,real(x), 'b'), hold on
stem(k,imag(x), 'r')
hold off
axis([0,8,-2,2])

sum(real(xr) - r.^k .* cos(+phi*k))
sum(imag(xr) - r.^k .* sin(+phi*k))
sum(real(xi) - r.^k .* cos(+phi*k))
sum(imag(xi) - r.^k .* sin(-phi*k))
sum(real(x) - r.^k .* sin(+phi*k) * 2)

b = [0 3/4*sqrt(3) 0];
a = [1 -3/4 9/16];
h = impz(b,a, length(k));
sum(h-x')

fvtool(b,a)

zplane(b,a)


