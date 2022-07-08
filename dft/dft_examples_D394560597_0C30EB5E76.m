% [Signals- & Systems](https://github.com/spatialaudio/signals-and-systems-exercises),
% [University of Rostock](https://www.uni-rostock.de/en/),
% [Institute of Communications Engineering](https://www.int.uni-rostock.de/),
% [Sascha Spors](https://orcid.org/0000-0001-7225-9992),
% [Frank Schultz](https://orcid.org/0000-0002-3010-0294),
% [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)

clear all;
clc;

N = 8;
k = [0:N-1].';
K = k.*k.'; % outer product
F = exp(+1j*2*pi/N * K);  %Fourier matrix

% matlab / python / engineering convention
DFT = F';
IDFT = 1/N*F;

%% Task 11.1 D394560597
phi = -pi/4
kappa = -phi/(2*pi/N) % positive is delay, negative is leading
 
% complex x[k]
xmu = [0; N*exp(1j*phi); 0; 0; 0; 0; 0; 0]
xk = IDFT*xmu

subplot(1,1,1)
plot(k, cos(2*pi/N*k), 'or:', 'markersize',8), hold on
plot(k, sin(2*pi/N*k), 'ob:', 'markersize',8)
plot(k, real(IDFT(:,2))*N, 'xr:', 'markersize',8)
plot(k, imag(IDFT(:,2))*N, 'xb:')

plot(k, real( IDFT(:,2) * xmu(2) ), 'or-')
plot(k, imag( IDFT(:,2) * xmu(2) ), 'ob-')

%using phi
plot(k, cos(2*pi/N*k + phi), 'xr', 'markersize',8)
plot(k, sin(2*pi/N*k + phi), 'xb', 'markersize',8)
% using kappa
plot(k, cos(2*pi/N*(k - kappa)), 'xr', 'markersize',8)
plot(k, sin(2*pi/N*(k - kappa)), 'xb', 'markersize',8)

hold off
grid on
xlabel('k')

% real x[k]
xmu = [0; N/2*exp(1j*phi); 0; 0; 0; 0; 0; N/2*exp(-1j*phi)]
xk = IDFT*xmu;
if allclose(imag(xk), zeros(N,1))
    xk = real(xk)
end


%% Task 11.2 0C30EB5E76
xk = exp(+1j*2*pi/N * 2.5 * k) % 2pi/(2N/5), sequence is thus 2N periodic,
% but not periodic in N

xmu = DFT*xk;
abs(xmu)
angle(xmu) / pi * 16

figure
subplot(2,1,1)
plot(k,abs(xmu), 'o-')
xlabel('\mu')
ylabel('|X|')
subplot(2,1,2)
plot(k,angle(xmu), 'o-')
xlabel('\mu')
ylabel('\angle X')


%##############################################################################
function flag = allclose(a, b)
% https://numpy.org/doc/stable/reference/generated/numpy.allclose.html
% numpy.allclose(a, b, rtol=1e-05, atol=1e-08, equal_nan=False)
% https://stackoverflow.com/questions/28975822/matlab-equivalent-for-numpy-allclose
rtol=1e-05;
atol=1e-08;
flag = all( abs(a(:)-b(:)) <= atol+rtol*abs(b(:)) );
end
