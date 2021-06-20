clear all
%close all
clc

% page 420, example 5 from
%@book{Strang2007,
%	address = {Wellesley},
%	author = {Gilbert Strang},
%	publisher = {Wellesley-Cambridge},
%	title = {Computational Science and Engineering},
%	year = {2007}}

% Laurent series |z|>0
N = 4;
h_tmp = factorial(0:N) .* (-1).^[0:N];
h_tmp = 1 ./ h_tmp;
h = zeros(1,2*N+1);
h(1:2:end) = h_tmp;  % non-zero FIR coefficients

%frequency vector:
NW = 2^10;
dW = 2*pi / NW;
W = 0:dW:2*pi-dW;

% DTFT frequency response of the FIR filter as H(z)
H = freqz(h,1,W);

% analytical function Ha(z) that belongs to above Laurent series
% cf. equation (5), page 420 in Strang2007
z = exp(1j*W);
Ha = exp(-z.^(-2));

sum(abs(Ha-H).^2) %check error

plot(W, abs(H), 'linewidth', 3), hold on
plot(W, abs(Ha), ':', 'linewidth', 3), hold off
xlabel('\Omega / rad')
ylabel('|H|')
grid on
xlim([0 2*pi])
ylim([0 3])
