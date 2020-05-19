clear all
close all
clc

full_analysis_flag = false

%% Bode Plot of LTI Systems
% Analysing Frequency Response of 1st/2nd Order LTI Systems
% make sure that you also check freqs(B, A, w), abs(), angle() for bode plot
% visualization manually

% Signals- & Systems
% University of Rostock, [Institute of Communications Engineering]
% (https://www.int.uni-rostock.de/)
% Prof. [Sascha Spors](https://orcid.org/0000-0001-7225-9992)
% [Frank Schultz](https://orcid.org/0000-0002-3010-0294)
% [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)
%%
% UE 5.1: lowpass 2nd order, PT2 element
sz = [0];
sp = [0, -3/4-1j, -3/4+1j];
k = 25/16;
sys = zpk(sz, sp, k)
txt = 'Conjugate Complex Pole, -3/4 \pm 1j'
lti_bode_plot(sys, txt, full_analysis_flag)
%%
% UE 5.2: bandpass 2nd order
sz = [0];
sp = [-1/10, -10];
k = 10;
sys = zpk(sz, sp, k)
txt = 'bandpass 2nd order'
lti_bode_plot(sys, txt, full_analysis_flag)
%%
% UE 5.3: resonsant highpass 2nd order
sz = [0, 10];
sp = [-1-1i, -1+1i]/sqrt(2);
k = 1;
sys = zpk(sz, sp, k)
txt = 'resonsant highpass 2nd order'
lti_bode_plot(sys, txt, full_analysis_flag)

%%
% UE 6.1: maximum phase
sz = [2];
sp = [-1/2];
k = 2;
sys = zpk(sz, sp, k)
txt = '1st order maximum phase shelving'
lti_bode_plot(sys, txt, full_analysis_flag)
%%
% UE 6.1: minimum phase
sz = [-2];
sp = [-1/2];
k = 2;
sys = zpk(sz, sp, k)
txt = '1st order minimum phase shelving'
lti_bode_plot(sys, txt, full_analysis_flag)
%%
% UE 6.1: allpass
sz = [+2];
sp = [-2];
k = 1;
sys = zpk(sz, sp, k)
txt = '1st order allpass'
lti_bode_plot(sys, txt, full_analysis_flag)



%%
function lti_bode_plot(sys, txt, full_analysis_flag)
    %disp('poles:')
    %pole(sys)
    %disp('zeros:')
    %zero(sys)
    w = logspace(-2,+2,2^10)';
    [mag, phase, ~] = bode(sys, w);
    mag_db = 20*log10(squeeze(mag));
    phase_deg = squeeze(phase);

    figure
    subplot(2,1,1)
    semilogx(w,mag_db,'linewidth',3)
    grid on
    xlabel('\omega / (rad/s)')
    ylabel('level 20lg |H| in dB')
    title(txt)
    subplot(2,1,2)
    semilogx(w,phase_deg,'linewidth',3)
    grid on
    xlabel('\omega / (rad/s)')
    ylabel('phase arg(H) in deg')
    if verLessThan('matlab','8.5')
        % -- Code to run in MATLAB R2015a and earlier here --
    else
        % -- Code to run in MATLAB R2015a and later here --
        % full LTI system analysis typically covers:
        if full_analysis_flag
            linearSystemAnalyzer({'impulse'; 'step'; 'bode';...
            'nyquist'; 'nichols'; 'pzmap'},sys)
        end
    end
end