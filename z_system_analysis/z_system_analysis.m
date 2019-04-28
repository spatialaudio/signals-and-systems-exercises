%[Signals- & Systems](https://github.com/spatialaudio/signals-and-systems-exercises),
%[University of Rostock](https://www.uni-rostock.de/en/),
%[Institute of Communications Engineering](https://www.int.uni-rostock.de/),
%Prof. [Sascha Spors](https://orcid.org/0000-0001-7225-9992),
%[Frank Schultz](https://orcid.org/0000-0002-3010-0294),
%Till Rettberg,
%[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)
clear all
close all
clc

%% coefficients of discrete-time LTI system
% see https://de.mathworks.com/help/signal/ref/freqz.html for convention of
% z-domain transfer function
%        b0 z^0 + b1 z^-1 + b2 z^-2 + ...
% H(z) = --------------------------------
%        a0 z^0 + a1 z^-1 + a2 z^-2 + ...
% with
% b = [b0, b1, b2], a = [a0, a1, a2]
% H(z) is often normalized such that a0=1

% check discrete-time systems that were previously discussed in the tutorials
if false % exercise 8.1
    b = [+1, -1, +2];
    a = [+1, -1/2, +1/4]; % note the sign reversal of a1, a2 
    % compared to block diagram
    fs = 1; % sampling frequency in Hz
    % unit sampling frequency fs = 1 is default in Matlab routines, 
    % such as in the fvtool(b,a)
    % in this case physical frequency f and W/(2*pi) are equal
    % for audio fs = 44.1, 48, 88.2, 96, 192 kHz are typical choices
elseif false % exercise 7.3: system H1
    % FIR filter with impulse response vector h would be handled as
    b = [1,1,-1,1/2]; % = h
    a = 1;
    fs = 1; % sampling frequency in Hz
elseif true % exercise 7.3: system H3
    % stable IIR filter 
    b = [2, 0, 1];
    a = [1, -1/2];
    fs = 1; % sampling frequency in Hz
elseif false % ODE RLC-example
    % digital filter design with so called bilinear transform of a 
    % continuous-time ODE (we will learn this in detail in the DSP course):
    % ODE RLC-example from 'solving_2nd_order_ode.pdf' /
    % 'frequency_response_2nd_order_ode.pdf' is
    % 16/25 y''(t) + 24/25 y'(t) + y(t) = DiracDelta(t), y'(t=0)=0, y(t=0)=0
    % for example sampled with
    fs = ceil((5/4)/(2*pi) * 100); % sampling frequency in Hz, note: omega0 = 5/4
    %note that we just round up to integer for nicer plotting
    [b, a] = bilinear([25/16], [1, 24/16, 25/16], fs); % same as:
    %[b, a] = bilinear([1], [16/25, 24/25, 1], fs);
    dBRange = [-80,10];
else % pass through system
    b = 1;
    a = 1;
    fs = 1;
end

disp(['b = ', num2str(b)])
disp(['a = ', num2str(a)])

sos = [b a]; %second order structur (sos) convention, see sos2tf(), then a0:=1
%this is very often used to cascade second order systems in series, cf. sosfilt()

%% helping stuff
N = 2^13;  % length of evaluated signals
k = [0:N-1]'; %sample index

% color handling like Python's matplotlib with Matlab's standard colors:
C0 = [0, 0.4470, 0.7410];
C1 = [0.8500, 0.3250, 0.0980];
C2 = [0.9290, 0.6940, 0.1250];

% open a figure with normalized size
figure('Units', 'normalized',...
    'Position', [0.25 0.1 0.55 0.75],...
    'Name', 'analysis of discrete-time LTI system')

%% discussion in discrete-time domain

% impulse respone
h = impz(b,a,N); % either
h1 = zeros(N,1); h1(1) = 1; h1 = filter(b, a, h1); % or filter the unit impulse

% step response
he = stepz(b,a,N); % either
he1 = cumsum(h);  % or cumulative sum
he2 = filter(b, a, ones(N,1));  % or filter the unit step signal
% for computer finite length signals are required, thus a rect of length N
% is used

%do the plot job
subplot(3,2,1)
stem(k, h, 'LineWidth', 2, 'Color', C0)
xlim([0, 64]) % hard coded for the above examples
xlabel('k')
ylabel('impulse response h[k]')
grid on

subplot(3,2,2)
stem(k, he, 'LineWidth', 2, 'Color', C0)
xlim([0, 64]) % hard coded for the above examples
xlabel('k')
ylabel('step response h\epsilon[k]')
grid on

%% discussion in (digital) frequency domain
H = fft(h); % get the DTFT frequency response as a DFT-approximation
% 
% the DFT is usually calculated by the Fast Fourier Transform (FFT)
% FFT subsumes algorithms that calculate DFT with very efficient computational
% effort
%
% this DFT handling works only if h includes the decaying tail of the impulse
% reponse with sufficient precision, thus it will NOT work for unstable systems
%
% instead of doing the DFT approach, one might use the function
% [H, W] = freqz(b,a,N, 'whole'); % or
% [H, f] = freqz(b,a,N,fs, 'whole');  %or
% [H] = freqz(b,a,f,fs);
% to evaluate frequency response of impulse response h / system with coeff b,a
% on the whole unit circle (i.e. evaluating one period of the DTFT)
%
% for freqz(b,a,f,fs) and for the DFT approach one will need a suitable
% frequency vector f with arbitrary high resolution set up by N:
df = fs/N;  % frequency step between DFT bins
f = [0:df:fs-df]'; % eigenfrequencies of DFT in Hz
% using this handling, one must not worry about odd/even N!
% this then includes DC up to fs-df, corresponding to W = [0...2pi)
%
% one should realize that DTFT -> DFT is a sampling process, having
% implications on the signal in time domain (aliasing, periodization) 

W = 2*pi*f/fs; % digital angular frequency in rad, note that discrete-time
% system does not require fs, since in the first instance 0<=W<2*pi does not 
% rely on it
% the way how we implemented this here shall help to grasp
% the link to physical frequency f

if exist('dBRange')==0 % if not defined above, we should do it here
    dBRange = [min(20*log10(abs(H))) max(20*log10(abs(H)))+eps]; %+eps needed
    %for the case that min/max yield exactly the same results
    %one might elaborate that this occurs for the example b=1, a=1
end

% do the plot job, we demonstrate frequency axis handling
% I. logarithmic x-axis f/Hz 
% II. linear x-axis W/2/pi = f/fs 
% III. linear x-axis W
% IV. linear x-axis W/pi (Matlab default in e.g. fvtool(b,a) and filter design)
% in the subplots below

% magnitude with new Matlab's plot objects handling for log/lin x-axis in one plot
subplot(3,2,3)
line(W/(2*pi), 20*log10(abs(H)),'Color',C1,'LineWidth',2)
ax1 = gca;
ax2 = axes('Position',ax1.Position,...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none');
line(f, 20*log10(abs(H)),'Parent',ax2,'Color',C0,'LineWidth',2)
ax1.XColor = C1;
ax2.XColor = C0;
ax1.YColor = C1;
ax2.YColor = C0;
ax1.XScale = 'lin';
ax2.XScale = 'log';
ax1.XLim = [0, 1];
ax2.XLim = [f(2), fs]; % we cannot plot f=0 in log, thus we start at f(2)=df
ax1.YLim = dBRange;
ax2.YLim = dBRange;
ax1.XLabel.String = '\Omega / (2\pi) = f / f_s';
ax2.XLabel.String = ['frequency in Hz   or   f / Hz   for f_s=', num2str(fs), ' Hz'];
ax1.YLabel.String = 'magnitude in dB   or   20 lg|H| / dB';
ax2.YLabel.String = 'magnitude in dB   or   20 lg|H| / dB';

% phase plot with Matlab's traditional plot calling functions
subplot(3,2,4)
plot(W, angle(H)*180/pi, 'LineWidth', 2, 'Color', C0)
xlim([0, 2*pi])
ylim([-180, +180])
xticks([0:pi/4:2*pi])
xticklabels({'0', '\pi/4', '\pi/2', '3/4\pi', '\pi', '5/4\pi', '3/2\pi', '7/4\pi', '2\pi'})
yticks([-180:45:+180])
xlabel('digital angular frequency in radian   or   \Omega / rad')
ylabel('phase in degree   or   \angle H / deg')
grid on

%% z-plane, i.e. pole/zero plot
subplot(3,2,5)
[hz, hp, ht] = zplane(b, a); %b,a must be row vectors!!!
hz.LineWidth = 2;
hz.MarkerSize = 9;
hz.Color = C0;
hp.LineWidth = 2;
hp.MarkerSize = 9;
hp.Color = C1;
grid on
xlabel('real part \Re(z)')
ylabel('imaginary part \Im(z)')
title('z-plane')
legend('zeros', 'poles')

%% group delay
gd = grpdelay(b, a, f, fs);
subplot(3,2,6)
plot(W/pi, gd/fs, 'LineWidth', 2, 'Color', C0)
xlabel('\Omega / \pi = f / (f_s/2)');
ylabel('groupdelay in seconds   or   \tau_{GD} / s')
grid on
xlim([0, 2])

%%
set(gcf,'Units', 'pixels') % restore default
