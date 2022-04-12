clear all
close all
clc

%% Bode Plot of LTI Systems
% Analysing Frequency Response of 2nd Order Ordinary Differential Equation
% (ODE) with Constant Coefficients
% make sure that you also check freqs(B, A, w), abs(), angle() for bode plot
% visualization manually

% Signals- & Systems
% University of Rostock, [Institute of Communications Engineering]
% (https://www.int.uni-rostock.de/)
% Prof. [Sascha Spors](https://orcid.org/0000-0001-7225-9992)
% [Frank Schultz](https://orcid.org/0000-0002-3010-0294)
% [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)

%% Unity Gain
sz = 0;
sp = 0;
k = +1;
txt = 'Unity Gain'
sys = zpk(sz, sp, k)  % Matlab function
lti_bode_plot(sys, txt);  % function defined at end of the script

%% Gain and Polarity
sz = 0;
sp = 0;
k = -10;
txt = '20 dB Gain with Inverted Polarity'
sys = zpk(sz, sp, k)
lti_bode_plot(sys, txt);

%% Poles / Zeros in Origin
sz = [0, 0];  % note: more zeros than poles is not a causal system!
sp = [0];
k = 1;
sys = zpk(sz, sp, k)
txt = [num2str(length(sz)), ' Zeros / ', num2str(length(sp)), ' Poles in Origin'];
txt1 = [': ', num2str((length(sz)-length(sp))*20), ' dB / decade']; 
lti_bode_plot(sys, [txt, txt1]);

sz = [0];  
sp = [0, 0, 0];  % more poles than zeros
k = 1;
sys = zpk(sz, sp, k)
txt = [num2str(length(sz)), ' Zeros / ', num2str(length(sp)), ' Poles in Origin'];
txt1 = [': ', num2str((length(sz)-length(sp))*20), ' dB / decade']; 
lti_bode_plot(sys, [txt, txt1]);
clear txt1

%% Single Real Pole
sz = [0];
sp = [0, -1];
k = 1;
sys = zpk(sz, sp, k)
txt = 'Single Real Pole, decreasing slope, -20 dB / decade'
lti_bode_plot(sys, txt)

%% Single Real Zero
sz = [0, -1];
sp = [0];
k = 1;
sys = zpk(sz, sp, k)
txt = 'Single Real Zero, increasing slope, + 20 dB / decade'
lti_bode_plot(sys, txt);

%% Conjugate Complex Zero
sz = [0, -3/4-1j, -3/4+1j];
sp = [0];
k = 16/25;
sys = zpk(sz, sp, k)
txt = 'Conjugate Complex Zero'
lti_bode_plot(sys, txt)

%% Conjugate Complex Pole

% this is the ODE RLC-example from 'solving_2nd_order_ode.pdf':
% 16/25 y''(t) + 24/25 y'(t) + y(t) = DiracDelta(t), y'(t=0)=0, y(t=0)=0

%B = [0, 0, 1];
%A = [16/25, 24/25, 1];
%sys = tf(B, A)
%txt = 'Conjugate Compley Pole, -3/4 \pm 1j'
%lti_bode_plot(sys, txt)
%equal to the zeros, poles and DC gain:
sz = [0];
sp = [0, -3/4-1j, -3/4+1j];
k = 25/16;
sys = zpk(sz, sp, k)
txt = 'Conjugate Complex Pole, -3/4 \pm 1j'
lti_bode_plot(sys, txt)
% w = 0 -> 0 dB:
% 20*log10((25/16)/(sqrt(9/16+(0+1)^2)*sqrt(9/16+(0-1)^2)))
% w = 1 -> -0.2169 dB:
% 20*log10((25/16)/(sqrt(9/16+(1+1)^2)*sqrt(9/16+(1-1)^2)))
% w = 10 -> -36.0865 dB:
% 20*log10((25/16)/(sqrt(9/16+(10+1)^2)*sqrt(9/16+(10-1)^2)))
% w = 100 -> -76.1232 dB:
% 20*log10((25/16)/(sqrt(9/16+(100+1)^2)*sqrt(9/16+(100-1)^2)))

if verLessThan('matlab','8.5')
    % -- Code to run in MATLAB R2015a and earlier here --
else
    % -- Code to run in MATLAB R2015a and later here --
    % full LTI system analysis typically covers:
%    linearSystemAnalyzer({'impulse'; 'step'; 'bode';...
%        'nyquist'; 'nichols'; 'pzmap'},sys)
end

%% Other Conjugate Complex Pole Examples
%other examples with poles approaching the imag axis (i.e. lower damping):
sz = [0];
sp = [0, -1/2-1j, -1/2+1j];
k = 5/4;
sys = zpk(sz, sp, k)
txt = 'Conjugate Complex Pole, -1/2 \pm 1j'
lti_bode_plot(sys, txt)

sz = [0];
sp = [0, -1/4-1j, -1/4+1j];
k = 17/16;
sys = zpk(sz, sp, k)
txt = 'Conjugate Complex Pole, -1/4 \pm 1j'
lti_bode_plot(sys, txt)

sz = [0];
sp = [0, -1/8-1j, -1/8+1j];
k = 65/64;
sys = zpk(sz, sp, k)
txt = 'Conjugate Complex Pole, -1/8 \pm 1j'
lti_bode_plot(sys, txt)

%% Example: Bandpass
%B = [0,100,0]
%A = [10, 101, 10]
%sys = tf(B, A)
%txt = 'Bandpass'
%lti_bode_plot(sys, txt)
%equal to
sz = [0];
sp = [-1/10, -10];
k = 10;
sys = zpk(sz, sp, k)
txt = 'Bandpass'
lti_bode_plot(sys, txt)

%%
function lti_bode_plot(sys, txt)
    %disp('poles:')
    %pole(sys)
    %disp('zeros:')
    %zero(sys)
    w = [1e-2:1e-2:1e2]';
    [mag, phase, wout] = bode(sys, w);
    mag_db = 20*log10(squeeze(mag));
    phase_deg = squeeze(phase);

    figure
    subplot(2,1,1)
    semilogx(w,mag_db,'linewidth',3)
    grid on
    xlabel('\omega / (rad/s)')
    ylabel('Magnitude: abs(H) / dB')
    xlim([1e-2 1e+2])
    ylim([-40 40])
    yticks([-40:10:40])
    title(txt)
    subplot(2,1,2)
    semilogx(w,phase_deg,'linewidth',3, 'color', 'red')
    grid on
    xlabel('\omega / (rad/s)')
    ylabel('Phase: arg(H) / deg')
    xlim([1e-2 1e+2])
    ylim([-225 +225])
    yticks([-180:45:+180])
       
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