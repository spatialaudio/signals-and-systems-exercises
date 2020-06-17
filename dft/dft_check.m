% Permuation matrix vs. DFT matrix
clear all;
clc;

N = 8;
W = exp(1j*2*pi/N);
k = 0:N-1;

%these two matrices are used for diagonalizing P
F = W.^(k' .* k) %DFT matrix
L = diag(W.^k) %eigenvalues as diag matrix

%for N=4
P = [0,1,0,0;
     0,0,1,0;
     0,0,0,1;
     1,0,0,0]
%for arbitrary N
v = zeros(1,N);
v(2) = 1;
P = toeplitz([v(1) fliplr(v(2:end))], v)

%check diagonalization of P
disp('F^-1*P*F - L = 0 ?')
F^-1*P*F - L

% Matlab's numeric solution, this is not ideal, since it is not sorted as our
% DFT above, but of course works as well
[V, D] = eig(P);
disp('V^-1*P*V - D = 0 ?')
V^-1*P*V - D


%plot all eigenvalues, i.e. DFT 'eigen'-frequencies on DTFT unit cirlce:
phi = 0:2*pi/64:2*pi;
plot(cos(phi), sin(phi), 'k'), hold on
plot(real(diag(L)), imag(diag(L)),'ok', 'markersize', 10)
hold off
axis square
xlabel('Re(\lambda)')
ylabel('Im(\lambda)')
grid

% S. 217 ML Strang:
c = zeros(1,N);
c = [1,-1,1,-1,1,-1,1,-1]
C = toeplitz([c(1) fliplr(c(2:end))], c)

[V,D] = eig(C)
F*C(1,:)'
