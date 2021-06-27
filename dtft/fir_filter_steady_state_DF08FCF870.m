clear all
%close all
clc

N = 8*12;
h = [1/2 0 -1/2];
k = 0:N-1;
Om0 = pi/4;
x = exp(1i*Om0*k);

y = conv(h,x,'full');
kf = 0:length(y)-1;

x = real(x);
y = real(y);

plot(k,x,'ko'), hold on
plot(kf,y*sqrt(2),'bo')
plot(k,x,'k')
plot(kf,y*sqrt(2),'b')
hold off
xlim([8*10, 8*10 + 8])
ylim([-1, 1])
legend('x[k]', 'y[k] * sqrt(2)')
grid


N = 2^16;
[pd, w] = phasedelay(h,1,N);
[gd, w] = grpdelay(h,1,N);
idx = find(w<Om0,1, 'last');
w(idx) / pi;

disp('phase delay in samples')
pd(idx)  % -1: output one sample earlier than input for single frequency
disp('group delay in samples')
gd(idx)  % +1: output is one sample later than input for frequency group

