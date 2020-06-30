%[Signals- & Systems](https://github.com/spatialaudio/signals-and-systems-exercises),
%[University of Rostock](https://www.uni-rostock.de/en/),
%[Institute of Communications Engineering](https://www.int.uni-rostock.de/),
%[Sascha Spors](https://orcid.org/0000-0001-7225-9992),
%[Frank Schultz](https://orcid.org/0000-0002-3010-0294),
%[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)

% discrete-time LTI system examples
% minimum phase system in series with specific allpass to obtain
% a) mixed phase system (some zeros out of unit circle)
% b) maximum phase system (all zeros outside the unit circle)
clear all
close all
clc

%% minimum phase
z = [1/2*exp(+1j*3*pi/4), 1/2*exp(-1j*3*pi/4), 1/2];
p = [3/4*1j, -3/4*1j, -1/3];
k = 2;
sos_min = zp2sos(z,p,k)
plot_dtlti(sos_min, 'minimum phase')

%% a) mixed phase
z = [2*exp(+1j*3*pi/4), 2*exp(-1j*3*pi/4), 1/2];
p = [3/4*1j, -3/4*1j, -1/3];
k = 1/2;
sos_mix = zp2sos(z,p,k)
plot_dtlti(sos_mix, 'mixed phase')

%allpass for mixed phase to get Hmin * Hall = Hmix
z = [2*exp(+1j*3*pi/4), 2*exp(-1j*3*pi/4), 0]; %max phase zeros from mixed phase
p = [1/2*exp(+1j*3*pi/4), 1/2*exp(-1j*3*pi/4), 0]; % we need two SOS for fvtool(sos), thus zero/pole in origin
k = 1/4;
sos_all_mix = zp2sos(z,p,k)
plot_dtlti(sos_all_mix, 'allpass for mixed phase')


%% b) maximum phase
z = [2*exp(+1j*3*pi/4), 2*exp(-1j*3*pi/4), 2];
p = [3/4*1j, -3/4*1j, -1/3];
k = -1/4;
sos_max = zp2sos(z,p,k)
plot_dtlti(sos_max, 'maximum phase')

%allpass for maximum phase to get Hmin * Hall = Hmax
z = [2*exp(+1j*3*pi/4), 2*exp(-1j*3*pi/4), 2]; %all zeros from max phase
p = [1/2*exp(+1j*3*pi/4), 1/2*exp(-1j*3*pi/4), 1/2];
k = -1/8;
sos_all_max = zp2sos(z,p,k)
plot_dtlti(sos_all_max, 'allpass for maximum phase')

%% Hmin * Hall in series
plot_dtlti([sos_min; sos_all_mix], 'check Hmix = Hmin * Hall')

%% Hmin * Hall in series
plot_dtlti([sos_min; sos_all_max], 'check Hmax = Hmin * Hall')




% compact function for full discrete-time LTI system analysis
% in time/frequency domain
function plot_dtlti(sos, fig_name)
if size(sos,1) == 1  % some routines below behave weird if only one biquad 
    sos(end+1,:) = [1,0,0,1,0,0] % thus put another thru biquad in series
end
figure('Name', fig_name)
dW = 2*pi/2^10;
W = 0:dW:2*pi-dW;
W(1) = eps; % W(1)=0 is not working sometimes
H = freqz(sos,W);
gd = grpdelay(sos, W);
subplot(3,2,1)
plot(W, db(abs(H)))
ylabel('dB')
title('Level')
grid on
subplot(3,2,3)
plot(W, unwrap(angle(H)))
ylabel('rad')
title('Phase')
grid on
subplot(3,2,5)
plot(W, gd)
xlabel('\Omega / rad')
ylabel('samples')
title('Group Delay')
grid on
subplot(3,2,2)
impz(sos)
grid on
subplot(3,2,4)
stepz(sos)
grid on
subplot(3,2,6)
[z,p,k] = sos2zp(sos);
zplane(z,p)
text(0.8,0.8, ['gain=',num2str(k)])
title('pole/zero/gain plot z-plane')
grid on
end
