%[Signals- & Systems](https://github.com/spatialaudio/signals-and-systems-exercises),
%[University of Rostock](https://www.uni-rostock.de/en/),
%[Institute of Communications Engineering](https://www.int.uni-rostock.de/),
%Prof. [Sascha Spors](https://orcid.org/0000-0001-7225-9992),
%[Frank Schultz](https://orcid.org/0000-0002-3010-0294),
%[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)
clear all
close all
clc
% Task 11.3 C8864C8D9F

%% Matlab routines in 'time' domain:
x = [-1; 2; 4];
h = [3; 1; 5];
cconv(x, h, 3) %cyclic in 3 samples
conv(x, h) %linear

x = [-1; 2; 4; 0; 0];
h = [3; 1; 5; 0; 0];
cconv(x, h, 5) %cyclic in 5 samples, linear conv due to long enough zeropadding

%% %3-cyclic conv with circular matrix multiplication 
x = [-1; 2; 4];
X = toeplitz(x, [x(1); flipud(x(2:end))]);  %put x as 1st col and fill cyclic
h = [3; 1; 5];
H = toeplitz(h, [h(1); flipud(h(2:end))]); %put h as 1st col and fill cyclic
Y = X*H;
ym = Y(:,1); % 1st col has the result
y = ifft(fft(x) .* fft(h));
allclose(y, ym)
%% %5-cyclic conv of zeropadded signals = linear conv
x = [-1; 2; 4; 0; 0];
X = toeplitz(x, [x(1); flipud(x(2:end))]); %put x as 1st col and fill cyclic
h = [3; 1; 5; 0; 0];
H = toeplitz(h, [h(1); flipud(h(2:end))]); %put h as 1st col and fill cyclic
Y = X*H;
ym = Y(:,1); % 1st col has the result
y = ifft(fft(x) .* fft(h));
allclose(y, ym)


%##############################################################################
function flag = allclose(a, b)
% https://numpy.org/doc/stable/reference/generated/numpy.allclose.html
% numpy.allclose(a, b, rtol=1e-05, atol=1e-08, equal_nan=False)
% https://stackoverflow.com/questions/28975822/matlab-equivalent-for-numpy-allclose
rtol=1e-05;
atol=1e-08;
flag = all( abs(a(:)-b(:)) <= atol+rtol*abs(b(:)) );
end
