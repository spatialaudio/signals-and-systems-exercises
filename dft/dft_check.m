% Permutation, Fourier and DFT matrix
% based on Gilbert Strang:
% 2016, Introduction to Linear Algebra, 5th, Wellesley Cambridge, chapter 8.3
% 2019, Linear Algebra and Learning from Data, Wellesley Cambridge, chapter IV.1
%
%[Signals- & Systems](https://github.com/spatialaudio/signals-and-systems-exercises),
%[University of Rostock](https://www.uni-rostock.de/en/),
%[Institute of Communications Engineering](https://www.int.uni-rostock.de/),
%[Sascha Spors](https://orcid.org/0000-0001-7225-9992),
%[Frank Schultz](https://orcid.org/0000-0002-3010-0294),
%[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)
clear all;
clc;
debug_flag = 0;

N = 4;
Sign = +1; %-1, +1 is engineering convention
W = exp(Sign * 1j*2*pi/N);  % plus sign!
k = [0:N-1]';
if debug_flag
    % create a circulant matrix from the row vector k.T
    C = toeplitz([k(1); flipud(k(2:end))], k)
end
K = k.*k'; % outer product

%these two matrices are used for diagonalizing P:
%F = W.^K %DFT matrix, this is inverse DFT with engineering sign convention
%L = diag(W.^k) %set up eigenvalues as diag matrix
% this is probably more efficient
F = exp(Sign * 1j*2*pi/N * K) 
L = diag(exp(Sign * 1j*2*pi/N*k))

if debug_flag
    disp('F*F''-N*eye(N)=0?')
    F*F'-N*eye(N)
    disp('F''*F-N*eye(N)=0?')
    F'*F-N*eye(N)
end

%permutation matrix, N=4
P = [0,1,0,0;
    0,0,1,0;
    0,0,0,1;
    1,0,0,0];
%permutation matrix, for arbitrary N built from toeplitz from row vector v.T
v = zeros(N,1);
v(2) = 1;
P = toeplitz([v(1); flipud(v(2:end))], v)

if debug_flag
    disp('singular values of P all 1:')
    svd(P) 
    %check diagonalization of P
    disp('F^-1 *P*F - L = 0 ?')
    (F^-1)*P*F - L
    disp('1/N*F'' *P*F - L = 0 ?')
    (1/N*F')*P*F - L
    %rewritten as
    disp('F*L* F^-1 - P = 0 ?')
    F*L*(F^-1) - P
    disp('F*L* 1/N*F''-P = 0 ?')
    F*L*(1/N*F')-P
end

if 0
    % Matlab's numeric solution, this is not ideal, since it is not sorted as
    % our DFT above, but of course works as well
    [V, D] = eig(P);
    disp('V^-1*P*V - D = 0 ?')
    V^-1*P*V - D
end

%plot all unit gain eigenvalues, i.e. DFT 'eigen'-frequencies on the unit circle:
Nc = 32;
phi = [0:Nc-1]*2*pi/Nc;
plot(cos(phi), sin(phi), 'k'), hold on
plot(real(diag(L)), imag(diag(L)),'ok', 'markersize', 10)
hold off
axis square
xlabel('Re(\lambda)')
ylabel('Im(\lambda)')
title(['\lambda^',num2str(N),'=1'])
grid

% adapted to our DFT/iDFT convention, x...time, X...frequency:
% X[mu] =     sum_{k= 0}^{N-1} x[k]  * exp(-j 2pi/N k mu)
% x[k]  = 1/N sum_{mu=0}^{N-1} X[mu] * exp(+j 2pi/N k mu)
Forward_DFT_Matrix = F'; % = conj(F) works as well, since F=F^T 
Inverse_DFT_Matrix = 1/N * F;
x = zeros(N,1); x(2)=1;  % set up dirac at k=1
X = Forward_DFT_Matrix * x;  % DFT, spectral analysis
xd = Inverse_DFT_Matrix * X; % inverse DFT, signal synthesis
disp('X - fft(x) = 0?')
X-fft(x)
disp('x - iDFT(DFT(x)) = 0?')
x-xd



%c = k;
%c = fft(exp(1j*0*2*pi/N*k))
c = fft(exp(1j*1*2*pi/N*k))
%c = fft(exp(1j*2*2*pi/N*k))
%c = fft(exp(1j*3*2*pi/N*k))
%create circulant matrix from the row vector c.T
C = toeplitz([c(1); flipud(c(2:end))], c);

if debug_flag
disp('F*c - ifft(c)*N = 0 ?')
F*c - ifft(c)*N
end
%eigenvalues of C are encoded in F*c
F*c % which is the inverse DFT, i.e. the time signal

%
if debug_flag
P2 = P*P;
L2 = diag(exp(Sign * 1j*2*pi/N*k * 2))
P3 = P*P*P;
L3 = diag(exp(Sign * 1j*2*pi/N*k * 3))
P4 = P*P*P*P;
L4 = diag(exp(Sign * 1j*2*pi/N*k * 4))
(F^-1)*P2*F - L2
(F^-1)*P3*F - L3
(F^-1)*P4*F - L4
end

disp('cyclic convolution:')

c = [-1,2,4]';
d = [3,1,5]';

C = toeplitz([c(1); flipud(c(2:end))], c)
D = toeplitz([d(1); flipud(d(2:end))], d)

N = 3; k = [0:N-1]'; K = k.*k'; F = exp(1j*2*pi/N * K);

A = C*D  % cyclic convolution as matrix mult
cconv(c,d,N)  % cyclic conv
1/N*F' * ( (F*c) .* (F*d) )  % fast cyclic convolution

F * A(1,:)'  % Fourier of cyclic conv result
(F*c) .* (F*d)  % element-wise mult of Fourier c and Fourier d



