%[Signals- & Systems](https://github.com/spatialaudio/signals-and-systems-exercises),
%[University of Rostock](https://www.uni-rostock.de/en/),
%[Institute of Communications Engineering](https://www.int.uni-rostock.de/),
%[Sascha Spors](https://orcid.org/0000-0001-7225-9992),
%[Frank Schultz](https://orcid.org/0000-0002-3010-0294),
%[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)
clear all
close all
clc

N = 9; % DFT length
k = (0:N-1).'; %time
mu = (0:N-1); %frequency
K = k*mu; %matrix = outer product

%% Fourier Matrix F is orthogonal
F = exp(1j*2*pi/N * K); % Fourier matrix
allclose(F'*F / N, eye(N))
allclose(F*F' / N, eye(N))
%% unitary Fourier matrix
% we can make F orthonormal by 1/sqrt(N) normalization
Fn = F/sqrt(N);
% Fn is a unitary matrix, for which
% A^-1 == A^H == A^*
% A^-1 A = A A^-1 = I
% holds, we check this numerically:
allclose(Fn^(-1) , Fn')  % inv vs complex-conjugate transpos
allclose(Fn^(-1) , conj(Fn)) % inv vs. complex-conjugate
allclose(Fn'*Fn, eye(N))  % orthonormal must yield identity,  left inverse
allclose(Fn*Fn', eye(N))  % orthonormal must yield identity,  right inverse

%% DC signal (all ones)
% DFT/IDFT with orthogonal Fourier matrix
x = ones(N,1);
X = conj(F) * x; % DFT
x_r = 1/N * F * X; % IDFT
allclose(x_r, x)

% DFT/IDFT with unitary Fourier matrix
x = ones(N,1);
X = conj(F) * x; % DFT
x_r = Fn * X;  % IDFT
allclose(x_r, x)

% fft()/ifft() vs. orthogonal matrix ops
allclose(fft(x) , conj(F) * x)
allclose(ifft(x), 1/N * F * x)

% fft()/ifft() vs. unitary matrix ops
allclose(fft(x)/sqrt(N) , conj(Fn) * x)
allclose(ifft(x)*sqrt(N), Fn * x)



%##############################################################################
function flag = allclose(a, b)
% https://numpy.org/doc/stable/reference/generated/numpy.allclose.html
% numpy.allclose(a, b, rtol=1e-05, atol=1e-08, equal_nan=False)
% https://stackoverflow.com/questions/28975822/matlab-equivalent-for-numpy-allclose
rtol=1e-05;
atol=1e-08;
flag = all( abs(a(:)-b(:)) <= atol+rtol*abs(b(:)) );
end
