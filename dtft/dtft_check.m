% Signals- & Systems
% University of Rostock, [Institute of Communications Engineering]
% (https://www.int.uni-rostock.de/)
% Prof. [Sascha Spors](https://orcid.org/0000-0001-7225-9992)
% [Frank Schultz](https://orcid.org/0000-0002-3010-0294)
% [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)

clear all
%close all
clc

%create the DTFT numerator by intentionally 10 simple zeros
%https://www.wolframalpha.com/input/?i=1+%2B+x%5E-2+%2B+17%2F4*x%5E-4+%2B+17%2F4*x%5E-6+%2B+x%5E-8+%2B+x%5E-10
z = sym('z');
r = 1/sqrt(2);
H =...
(z-r*exp(+1i*1/4*pi))*...
(z-r*exp(-1i*1/4*pi))*...    
(z-1/r*exp(+1i*1/4*pi))*...
(z-1/r*exp(-1i*1/4*pi))*...
(z-r*exp(+1i*3/4*pi))*...
(z-r*exp(-1i*3/4*pi))*...    
(z-1/r*exp(+1i*3/4*pi))*...
(z-1/r*exp(-1i*3/4*pi))*...
(z-1i)*(z+1i);
expand(H * z^-10)
expand(H)
factor(H)
factor(1/z^10) % 10 poles in origin

%%

% corresponding FIR
h = [1 0 1 0 17/4 0 17/4 0 1 0 1] % nice coeff, nice zeros
% step response
he = cumsum(h)
disp('sum(h) is DC gain:')
sum(h)
disp('DC gain in dB')
20*log10(sum(h))

% check the DTFT, we use the FFT to evaluate the unit circle
N = 2^9;
H = fft(h,N);
Om = 2*pi/N*[0:N-1];

%%
dB = [-25 25];
subplot(6,2,1)
plot(Om,20*log10(abs(H)))
xlabel('\Omega / (rad/s)')
ylabel('dB')
axis([0 2*pi dB])

subplot(6,2,2)
plot(Om,unwrap(angle(H))), hold on
plot([pi/2 pi/2], [0 -30])
plot([2*pi-pi/2 2*pi-pi/2], [0 -30])
hold off
xlabel('\Omega / (rad/s)')
ylabel('rad')

cos0 = h(1) * exp(-1j*Om*0) + h(11) * exp(-1j*Om*10);
cos1 = h(2) * exp(-1j*Om*1) + h(10) * exp(-1j*Om*9);
cos2 = h(3) * exp(-1j*Om*2) + h(9) * exp(-1j*Om*8);
cos3 = h(4) * exp(-1j*Om*3) + h(8) * exp(-1j*Om*7);
cos4 = h(5) * exp(-1j*Om*4) + h(7) * exp(-1j*Om*6);
cos5 = h(6) * exp(-1j*Om*5);

if 1
    cos0 = cos0.*exp(1i*Om*5);
    cos1 = cos1.*exp(1i*Om*5);
    cos2 = cos2.*exp(1i*Om*5);
    cos3 = cos3.*exp(1i*Om*5);
    cos4 = cos4.*exp(1i*Om*5);
    cos5 = cos5.*exp(1i*Om*5);
end

subplot(6,2,3)
plot(Om,real(cos0)), hold on
plot(Om,imag(cos0))
plot(Om, abs(cos0))
hold off

subplot(6,2,4)
plot(Om,real(cos1)), hold on
plot(Om,imag(cos1))
plot(Om, abs(cos1))
hold off

subplot(6,2,5)
plot(Om,real(cos2)), hold on
plot(Om,imag(cos2))
plot(Om, abs(cos2))
hold off

subplot(6,2,6)
plot(Om,real(cos3)), hold on
plot(Om,imag(cos3))
plot(Om, abs(cos3))
hold off

subplot(6,2,7)
plot(Om,real(cos4)), hold on
plot(Om,imag(cos4))
plot(Om, abs(cos4))
hold off


subplot(6,2,8)
plot(Om,real(cos5)), hold on
plot(Om,imag(cos5))
plot(Om, abs(cos5))
hold off


Hman = real(cos0) + real(cos2) + real(cos4) + ...
    1i*imag(cos0) + 1i*imag(cos2) + 1i*imag(cos4);
%Hman = cos0 + cos2 + cos4;
Hmanabs = abs(Hman);
HmandB = 20*log10(Hmanabs);

disp('check that imag is zero')
min(imag(Hman))
max(imag(Hman))

subplot(6,2,9)
plot(Om,real(Hman)), hold on
plot(Om,Hmanabs)
plot([pi/2 pi/2], [-4 4])
plot([2*pi-pi/2 2*pi-pi/2], [-4 4])
hold off
xlabel('\Omega / (rad/s)')

subplot(6,2,10)
plot(Om,HmandB, 'linewidth',3), hold on
plot(Om,20*log10(abs(H)))
plot([pi/2 pi/2], dB)
plot([2*pi-pi/2 2*pi-pi/2], dB)
hold off
axis([0 2*pi dB])
xlabel('\Omega / (rad/s)')

subplot(6,2,11)
zplane(h,1)

subplot(6,2,12)
stem([0:length(h)-1], h), hold on
stem([0:length(he)-1], he), hold off
xlabel('k')


