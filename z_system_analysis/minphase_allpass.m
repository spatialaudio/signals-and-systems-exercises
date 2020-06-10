%[Signals- & Systems](https://github.com/spatialaudio/signals-and-systems-exercises),
%[University of Rostock](https://www.uni-rostock.de/en/),
%[Institute of Communications Engineering](https://www.int.uni-rostock.de/),
%Prof. [Sascha Spors](https://orcid.org/0000-0001-7225-9992),
%[Frank Schultz](https://orcid.org/0000-0002-3010-0294),
%Till Rettberg,
%[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)
% UE 8.1
clear all
close all
clc

% color handling like Python's matplotlib:
C0 = [0, 0.4470, 0.7410];
C1 = [0.8500, 0.3250, 0.0980];
C2 = [0.9290, 0.6940, 0.1250];

%% system
b = [+1 -1 +2];
a = [+1 -1/2 +1/4];
[z,p,k] = tf2zpk(b,a)
subplot(1,4,1)
zplane(b,a)
xlim([-1.5 +1.5])
ylim([-1.5 +1.5])
grid on
title('H(z), k=1')
subplot(1,4,4)
zplane(b,a)
hold on

%% minphase system
bMin = [+1 -1/2 +1/2]; % inverted zeros of H(z), 1./z == zMin
aMin = [+1 -1/2 +1/4]; % poles of H(z)
[zMin,pMin,kMin] = tf2zpk(bMin,aMin)
subplot(1,4,2)
zplane(bMin,aMin)
xlim([-1.5 +1.5])
ylim([-1.5 +1.5])
grid on
title('H_{Min}(z), k=1')
subplot(1,4,4)
zplane(bMin,aMin)
hold on

%% allpass system
bAll = [+1 -1 +2]; % zeros of H(z)
aAll = [+1 -1/2 +1/2]; %zeros of minphase
[zAll,pAll,kAll] = tf2zpk(bAll,aAll)
subplot(1,4,3)
zplane(bAll,aAll), hold on
plot([0 real(zAll(1))],[0 imag(zAll(1))], 'color', C1)
plot([0 real(zAll(2))],[0 imag(zAll(2))], 'color', C1)
xlim([-1.5 +1.5])
ylim([-1.5 +1.5])
grid on
title('H_{All}(z), k=1')
subplot(1,4,4)
zplane(bAll,aAll)
hold off
xlim([-1.5 +1.5])
ylim([-1.5 +1.5])
grid on
title('H_{All}(z) \cdot H_{Min}(z) = H(z)k, k=1')
