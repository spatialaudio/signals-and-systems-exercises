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
k = (0:N-1)'; %time
mu = (0:N-1); %frequency
K = k*mu; %matrix = outer product
F = exp(1j*2*pi/N * K); % Fourier matrix

% F is orthogonal
norm(F'*F - N*eye(N))

% we can make F orthonormal by a simple 1/sqrt(N) normalization
F_unitary = F/sqrt(N);
% F_unitary is a unitary matrix, for which
% A^-1 == A^H == A^*
% A^H A = I
% holds, we check this numerically:
norm(F_unitary^(-1) - F_unitary')  % inv vs hermitian
norm(F_unitary^(-1) - conj(F_unitary)) % inv vs. complex-conjugate
norm(F_unitary'*F_unitary - eye(N))  % orthonormal must yield identity

% DFT/IDFT for DC signal (all ones)
x = ones(N,1);
X = F' * x; % DFT
x_r = 1/N * F * X; % IDFT
norm(x_r-x)

% DFT/IDFT with unitary matrix
x = ones(N,1);
X = F_unitary' * x; % DFT
x_r = F_unitary * X;  % IDFT
norm(x_r-x)

% fft()/ifft() vs. unitary matrix ops
norm((fft(x) / sqrt(N)) - (F_unitary' * x))
norm((ifft(x) * sqrt(N)) - (F_unitary * x))
