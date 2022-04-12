%[Signals- & Systems](https://github.com/spatialaudio/signals-and-systems-exercises),
%[University of Rostock](https://www.uni-rostock.de/en/),
%[Institute of Communications Engineering](https://www.int.uni-rostock.de/),
%[Sascha Spors](https://orcid.org/0000-0001-7225-9992),
%[Frank Schultz](https://orcid.org/0000-0002-3010-0294),
%[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)

% exercise 12.1 #6337B75DF2

% Frequency Response of FIR Filter
% a) analytical DTFT
% b) DFT
% c) zeropadded DFT
% d) DFT -> DTFT interpolation
clear all;
%close all;
clc;

% finite impulse response starting at k=0, length = 11, FIR type I
h = [1 0 1 0 17/4 0 17/4 0 1 0 1];

%% a) DTFT 
Ndtft = 2^10;
dW = 2*pi/Ndtft;
Wdtft = 0:dW:2*pi-dW;
Hdtft = (2*cos(5*Wdtft) + 2*cos(3*Wdtft) + 17/2*cos(Wdtft)) .* exp(-1j*5*Wdtft);

%% b) DFT
Ndft = length(h);
dW = 2*pi/Ndft;
Wdft = 0:dW:2*pi-dW;
Hdft = fft(h);

%% c) zeropadded DFT
hz = zeros(1,16);
hz(1:Ndft) = h;
Ndftz = length(hz);
dW = 2*pi/Ndftz;
Wdftz = 0:dW:2*pi-dW;
Hdftz = fft(hz);
clear dW; 

%% d) DFT -> DTFT interpolation
Hint = interpolate_dft2dtft(Hdft, Wdtft);

%% plot
subplot(2,1,1)
stem(Wdft, abs(Hdft), 'linewidth', 2), hold on
plot(Wdtft, abs(Hdtft), 'linewidth', 2)
stem(Wdftz, abs(Hdftz), ':', 'linewidth', 2)
plot(Wdtft, abs(Hint), ':', 'linewidth', 2)
hold off
legend('N=11 DFT', 'DTFT', 'zeropadded N=16 DFT', 'DFT -> DTFT interpolation')
xlabel('\Omega / rad')
ylabel('magnitude')
title('Frequency Response of Linear Phase Type I FIR Filter with 11 Coefficients')
grid on

subplot(2,1,2)
stem(Wdft, unwrap(angle(Hdft)), 'linewidth', 2), hold on
plot(Wdtft, unwrap(angle(Hdtft)), 'linewidth', 2)
stem(Wdftz, unwrap(angle(Hdftz)), ':', 'linewidth', 2)
plot(Wdtft, unwrap(angle(Hint)), ':', 'linewidth', 2)
hold off
xlabel('\Omega / rad')
ylabel('unwrapped phase in rad')
grid on

%%
function Xint = interpolate_dft2dtft(X, W)
N = length(X);
tmp_2piN = 2*pi/N;
tmp_N2 = (N-1)/2;
Xint = W*0;

for Widx = 1:length(W)
    Xint(Widx) = 0;
    for mu=0:N-1
        W_tmp = W(Widx) - tmp_2piN*mu;
        Xint(Widx) = Xint(Widx) + X(mu+1)*diric(W_tmp, N)*exp(-1j*W_tmp*tmp_N2);
    end
end

end

% ## Copyright
% 
% This tutorial is provided as Open Educational Resource (OER), to be found at
% https://github.com/spatialaudio/signals-and-systems-exercises
% accompanying the OER lecture
% https://github.com/spatialaudio/signals-and-systems-lecture.
% Both are licensed under a) the Creative Commons Attribution 4.0 International
% License for text and graphics and b) the MIT License for source code.
% Please attribute material from the tutorial as *Frank Schultz,
% Continuous- and Discrete-Time Signals and Systems - A Tutorial Featuring
% Computational Examples, University of Rostock* with
% ``github URL, commit number and/or version tag, year, (file name and/or content)``.