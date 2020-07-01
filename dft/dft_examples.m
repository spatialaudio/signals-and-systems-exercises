clear all;
clc;
debug_flag = 1;

N = 8;
k = [0:N-1]';
K = k.*k'; % outer product
F = exp(+1j*2*pi/N * K);  %Fourier matrix

%matlab / python / engineering convention
DFT = F';
IDFT = 1/N*F;

% Task 11.1 D394560597
phi = -pi/4
kappa = -phi/(2*pi/N) % positive is delay, negative leading
 % ' operator is not good here, .' would be needed, so make it more obvious:

 
% complex x[k]
xmu = transpose([0 N*exp(1j*phi) 0 0 0 0 0 0])
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
xmu = transpose([0 N/2*exp(1j*phi) 0 0 0 0 0 N/2*exp(-1j*phi)])
xk = IDFT*xmu




%%
% Task 11.2 0C30EB5E76
xk = exp(+1j*2*pi/N * 2.5 * k) % 2pi/(2N/5), thus 2N periodic

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
