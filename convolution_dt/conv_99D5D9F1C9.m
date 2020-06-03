clear all
%close all
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




