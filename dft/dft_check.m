% Permutation, Fourier and DFT matrix
% based on Gilbert Strangs:
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
W = exp(+1j* Sign *2*pi/N);  % +1j!
k = [0:N-1].';

if debug_flag
    % create a circulant matrix from the row vector k.'
    C = toeplitz([k(1); flipud(k(2:end))], k)
end
K = k * k.'; % outer product

%these two matrices are used for diagonalizing permutation matrix P:
%F = W.^K %DFT matrix, this is inverse DFT when using engineering sign convention
%L = diag(W.^k) %set up eigenvalues as diag matrix
% this is probably more efficient to compute
F = exp(+1j* Sign *2*pi/N * K) 
L = diag(exp(+1j* Sign *2*pi/N*k))

if debug_flag
    disp('F*F''/N == eye(N)')
    allclose(F*F' / N, eye(N))
    disp('F''*F/N == eye(N)')
    allclose(F'*F/N, eye(N))
end

%a permutation matrix for N=4
P = [0,1,0,0;
    0,0,1,0;
    0,0,0,1;
    1,0,0,0];
%this permutation but for arbitrary N, built from toeplitz using row vector v.T
v = zeros(N,1);
v(2) = 1;
P = toeplitz([v(1); flipud(v(2:end))], v);

if debug_flag
    disp('singular values of P all 1:')
    svd(P) 
    %check diagonalization of P
    disp('F^-1*P*F == L')
    allclose((F^-1)*P*F , L)
    disp('1/N*F''*P*F, L')
    allclose((1/N*F')*P*F, L)
    %rewritten as
    disp('F*L*F^-1==P')
    allclose(F*L*(F^-1), P)
    disp('F*L*1/N*F''==P')
    allclose(F*L*(1/N*F'),P)
end

if 0
    % Matlab's numeric solution, this is not ideal, since it is not sorted as
    % our DFT above, but of course works as well
    [V, D] = eig(P);
    disp('V^-1*P*V == D')
    allclose(V^-1*P*V , D)
end

%plot all unit gain eigenvalues, i.e. DFT 'eigen'-frequencies on the unit circle:
Nc = 2^8;
phi = [0:Nc-1]*2*pi/Nc;
plot(cos(phi), sin(phi), 'k'), hold on
plot(real(diag(L)), imag(diag(L)),'ok', 'markersize', 10)
hold off
axis square
xlabel('Re(\lambda)')
ylabel('Im(\lambda)')
title(['\lambda^',num2str(N),'=1'])
grid

disp('our convention, sign=+1')
% adapted to our DFT/iDFT convention, x...time, X...frequency:
% X[mu] =     sum_{k= 0}^{N-1} x[k]  * exp(-j 2pi/N k mu)
% x[k]  = 1/N sum_{mu=0}^{N-1} X[mu] * exp(+j 2pi/N k mu)
Forward_DFT_Matrix = conj(F); % = F' works as well, since F=F^T 
Inverse_DFT_Matrix = 1/N * F;
x = zeros(N,1); x(2)=1;  % set up Dirac at k=1
X = Forward_DFT_Matrix * x;  % DFT, signal analysis
xd = Inverse_DFT_Matrix * X; % inverse DFT, signal synthesis
if Sign==+1
    disp('X == fft(x)')
    allclose(X, fft(x))
end
disp('x == iDFT(DFT(x))')
allclose(x, xd)


%c = exp(1j*0*2*pi/N*k);
c = exp(1j*1*2*pi/N*k);
%c = exp(1j*2*2*pi/N*k);
%c = exp(1j*3*2*pi/N*k);
%create circulant matrix from the row vector c.T
C = toeplitz([c(1); flipud(c(2:end))], c);

if debug_flag
    if Sign==+1
        disp('F*c == fft(c)')
        allclose(F'*c , fft(c))
    elseif Sign==-1
        disp('F*c == ifft(c)*N')
        allclose(F'*c , ifft(c)*N)
    end
end
%eigenvalues of C are encoded in F'*c
if Sign==+1
    allclose(F'*c, ... % which is the DFT, i.e. a spectrum
    fft(c))
elseif Sign==-1
    allclose(F'*c, ... % which is the DFT, i.e. a spectrum
    N*ifft(c))
end

%
if debug_flag
    P2 = P*P;
    L2 = diag(exp(+1j* Sign *2*pi/N*k * 2))
    P3 = P*P*P;
    L3 = diag(exp(+1j * Sign *2*pi/N*k * 3))
    P4 = P*P*P*P;
    L4 = diag(exp(+1j * Sign *2*pi/N*k * 4))
    allclose((F^-1)*P2*F, L2)
    allclose((F^-1)*P3*F, L3)
    allclose((F^-1)*P4*F, L4)
end


disp('cyclic convolution:')

x = [-1; 2; 4];
h = [3; 1; 5];

X = toeplitz(x, [x(1); flipud(x(2:end))])
H = toeplitz(h, [h(1); flipud(h(2:end))])
Y = X*H  % cyclic convolution as matrix mult

N = 3;
k = [0:N-1]';
K = k*k.';
F = exp(+1j*2*pi/N * K);

allclose(cconv(x,h,N),...  % cyclic conv
1/N*F * ( (F'*x) .* (F'*h) ))  % vs fast cyclic convolution

allclose(F' * Y(:,1),...  % DFT of cyclic conv result
(F'*x) .* (F'*h))  % element-wise mult of DFT coeff


%##############################################################################
function flag = allclose(a, b)
% https://numpy.org/doc/stable/reference/generated/numpy.allclose.html
% numpy.allclose(a, b, rtol=1e-05, atol=1e-08, equal_nan=False)
% https://stackoverflow.com/questions/28975822/matlab-equivalent-for-numpy-allclose
rtol=1e-05;
atol=1e-08;
flag = all( abs(a(:)-b(:)) <= atol+rtol*abs(b(:)) );
end
