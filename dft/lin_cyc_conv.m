%[Signals- & Systems](https://github.com/spatialaudio/signals-and-systems-exercises),
%[University of Rostock](https://www.uni-rostock.de/en/),
%[Institute of Communications Engineering](https://www.int.uni-rostock.de/),
%Prof. [Sascha Spors](https://orcid.org/0000-0001-7225-9992),
%[Frank Schultz](https://orcid.org/0000-0002-3010-0294),
%[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)
clear all
close all
clc

%% Matlab routines in 'time' domain:
x = [-1 2 4]';
h = [3 1 5]';
conv(x,h) %linear
cconv(x,h,3) %cyclic in 3 samples

x = [-1 2 4 0 0]';
h = [3 1 5 0 0]';
cconv(x,h,5) %cyclic in 5 samples, linear conv due to zeropadding

%% %3-cyclic conv
x = [-1 2 4]';
X = toeplitz([x(1); flipud(x(2:end))], x);  %put x as 1st row and fill cyclic
h = [3 1 5]';
H = toeplitz([h(1); flipud(h(2:end))], h); %put h as 1st row and fill cyclic
Y = X*H;
ym = Y(1,:)' % 1st row has the result
y = ifft(fft(x) .* fft(h));
y-ym
%% %5-cyclic conv of zeropadded signals = linear conv
x = [-1 2 4 0 0]';
X = toeplitz([x(1); flipud(x(2:end))], x); %put x as 1st row and fill cyclic
h = [3 1 5 0 0]';
H = toeplitz([h(1); flipud(h(2:end))], h); %put h as 1st row and fill cyclic
Y = X*H;
ym = Y(1,:)' % 1st row has the result
y = ifft(fft(x) .* fft(h));
y-ym
