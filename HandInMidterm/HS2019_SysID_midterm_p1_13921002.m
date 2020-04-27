function [p1_R,p1_omega,p1_a,p1_var] = HS2019_SysID_midterm_p1_13921002()

%% Solution for Problem 1
%% Output format specification
% p1_R must be a 1xT vector
% p1_omega must be a 1xM vector
% p1_a must be a 1xM vector
% p1_var must be a scalar
%% Generate data

% Extract Legi from Filename
name=mfilename;
LegiNumber= name(end-7:end);

p1_U = HS2019_SysID_midterm_p1_GenerateData(LegiNumber);

%% General instructions for solution

% Change the filename of this function, both in the function definition
% above and in the filename in the folder

% Use the variable p1_U to solve the problem. 

% Modify your code in the next sections, and return the variables
% requested.

% If you skip one part of the problem, return the empty vectors as already
% provided in the code

%% Task 1: Calculation of Autocorrelation
fprintf('\n')
fprintf('----------------------------------------------------------------------\n')
fprintf('----------------------------------------------------------------------\n')
fprintf('\n')
fprintf('Problem 1')
fprintf('\n')
fprintf('\n')
fprintf('----------------------------------------------------------------------\n')
fprintf('\n')
fprintf('Task 1: Calculation of periodic autocorrelation via Window 1')
fprintf('\n')
fprintf('\n')

Ts = 0.5;
K = length(p1_U);
t = (Ts * 0:(K-1))';
exp1 = p1_U(1,:);

% Use autocorrelation. I should get high correlation once the lag time matches the period.
ac_tot = zeros(K,1);
lags_tot = zeros(K,1);
[ac_tot, lags_tot] = R_per(exp1, K);

% Plot autocorr of whole signal
figure(1)
clf
hold on
title('Task 1: Autocorrelation of the whole 1st-window signal', 'FontSize', 12)
stem(lags_tot, ac_tot)
xlim([0 K])
ylim([-80 150])
xlabel('\tau', 'FontSize',20)
ylabel('R(\tau)')

% Select period T from graphfprintf('From Figure 1 one can see that the autocorrelation plot shows\nenergy peaks every 280 values, so T = 280*Ts = 140s.')
M = 280;
T = M*Ts;

% Calculate the Periodic Autocorrelation of only one period.
exp1_per = exp1(1:M);
[ac_per, lags_per] = R_per(exp1_per, M);

%plot one period of the Autocorrelation Function
figure(2)
clf
hold on
title('Task 1: One period of the periodic autocorrelation function', 'FontSize', 12)
stem(lags_per(1.5*M+1:2.5*M), ac_per(1.5*M+1:2.5*M))
xlim([-M/2+1 M/2])
ylim([-70 150])
xlabel('\tau', 'FontSize', 20)
ylabel('R(\tau)','FontSize', 20)

%Select just the values between lags (-T/2+1) and (T/2)
ac_per_sel = ac_per(421:700);

% Explanation
fprintf('Autocorrelation peaks will appear after each period, so I plot\nthe autocorrelation of the whole signal and analyze the graph to\ndetect those peaks [Figure 1].\n')
fprintf('The peaks occur every 280 values, which leads to T=140s (280 samples)\ndue to Ts=0.5s.\n')
fprintf('\n')
fprintf('In order to compute the periodic autocorrelation vector p1_R,\nI use the formula presented during lecture in slide 2.19.\n')
fprintf('\n')
fprintf('One period of the autocorrelation function is displayed in Figure 2.')
fprintf('\n')
%fprintf('???INCLUDE??? Lags -139 to 140 correspond to the values stored by Matlab in\n???INCLUDE??? cells 421 to 700.\n')

p1_R = [transpose(ac_per_sel)];


%% Task 2: Estimation of signal components
fprintf('\n')

% Calculate how many points you can throw away
throw_points = mod(K,T);

% Define a new df p1_U_clean, without the first 112 points of each row
p1_U_clean = p1_U(:,throw_points+1:K);

% Calculate FFT of each row independently, since there could be a
% phase-shift.
WINDOWS_MATRIX = fft(p1_U_clean, length(p1_U_clean), 2);

WINDOWS_SUM = zeros(1, length(p1_U_clean));
for i = 1:5
    for n = 1:length(p1_U_clean)
        WINDOWS_SUM(1,n) = WINDOWS_SUM(1,n) + abs(WINDOWS_MATRIX(i,n));
    end
end    
WINDOWS_AVG = WINDOWS_SUM ./ 5;

%Calculate omegas
omega = (2*pi/(Ts*length(p1_U_clean)))*transpose([0:length(p1_U_clean)-1]);
idx = find(omega > 0 & omega < pi/Ts);

% Plot FFT single experiments
figure(3)
clf
hold on
title('Task 2: FFT of individual windows','FontSize', 12)
plot(omega(idx),(abs(WINDOWS_MATRIX(1,idx)) / length(p1_U_clean)) * 2,'blue')
plot(omega(idx),(abs(WINDOWS_MATRIX(2,idx)) / length(p1_U_clean)) * 2,'red')
plot(omega(idx),(abs(WINDOWS_MATRIX(3,idx)) / length(p1_U_clean)) * 2,'green')
plot(omega(idx),(abs(WINDOWS_MATRIX(4,idx)) / length(p1_U_clean)) * 2,'cyan')
plot(omega(idx),(abs(WINDOWS_MATRIX(5,idx)) / length(p1_U_clean)) * 2,'magenta')
legend('Window 1','Window 2','Window 3','Window 4','Window 5')
xlabel('\omega','FontSize', 20);
ylabel('|U(e^{jw})|','FontSize', 20);
xlim([0 pi/Ts])
ylim([0 10])

% Plot FFT average
figure(4)
clf
hold on
title('Task 2: U of averaged period','FontSize', 12)
stem(omega(idx),(abs(WINDOWS_AVG(idx)) / length(p1_U_clean)) * 2)
xlabel('\omega','FontSize', 20);
ylabel('|U(e^{jw})|','FontSize', 20);
xlim([0 pi/Ts])
ylim([0 10])

% Find indices with magnitudes bigger than 130.
idx = find(abs(WINDOWS_AVG(1,:))>650);
freqs = 2*pi/(Ts*length(p1_U_clean)) * (idx - ones(1,length(idx)));
% Only take the positive frequencies
sol_freqs = freqs(:,1:6);

% Find magnitudes that correspond to those indices and take the absolute and
% divide by "length(p1_U_clean)" to get the indices in the time domain.
magnitudes_abs = (abs(WINDOWS_AVG(idx)) / length(p1_U_clean)) * 2;

% Only take the magnitudes that correspond to the positive frequencies
sol_magnitudes_abs = magnitudes_abs(:,1:6);

% Explain
fprintf('----------------------------------------------------------------------\n')
fprintf('\n')
fprintf('Task 2: Estimation of p1_omega and p1_a')
fprintf('\n')
fprintf('\n')
fprintf('I do the FFT of each window independently without regarding the first\n112 samples, so that I get an integer number of periods, in this case\n6 full periods. ')
fprintf('That way I have a high frequency resolution.\n')
fprintf('Then I calculate the average of all 5 rows to get rid of as much\nnoise as possible.\n')
fprintf('\n')
fprintf('Figure 3 shows the absolute value of the FFT for each individual\nwindow over the positive frequencies.\n')
fprintf('Figure 4 shows the absolute value of the averaged FFT.\n')
fprintf('\n')
fprintf('From the averaged FFT plot [Figure 4], I got the frequencies\nand magnitudes over a threshold to get rid of the frequencies due to\nnoise.\n')
fprintf('\n')
fprintf('The formula I used to get magnitudes of the coefficients p1_a in the\ntime domain:\np1_a[i] = 2*abs(X[i])/N, where\nX[k]: Fourier coefficient in frequency domain\nN: total number of points of the signal')

p1_omega = [sol_freqs];
p1_a = [sol_magnitudes_abs];


%% Task 3: Estimation of noise variance
fprintf('\n')
fprintf('----------------------------------------------------------------------\n')
fprintf('\n')
fprintf('Task 3: Estimation of noise variance')
fprintf('\n')
fprintf('\n')

% Create u_real signal using the frequencies and magnitudes you found out in Task2 
per_true_abs = zeros(1,length(p1_U_clean));

for i = 1:6
    for n = 1:length(p1_U_clean)
        per_true_abs(1,n) = per_true_abs(1,n) + sol_magnitudes_abs(1,i) * cos(sol_freqs(1,i)*n);
    end    
end

t = 0:Ts:(length(p1_U_clean)/2-Ts);

% Plot u_real
figure(5)
clf
hold on
plot(t,per_true_abs)
plot(t, p1_U_clean(1,:))
title('Task 3: Reconstructed signal (no noise) vs. original noisy signal (Window 1)', 'FontSize', 12)
legend('Reconstructed singal with no noise','Original signal (Window 1)')
xlabel('t [sec]','FontSize',20);
ylabel('signal','FontSize',20);
xlim([0 length(p1_U_clean)/2])
ylim([-30 40])

% Set frequencies over threshold to zero
NOISE = WINDOWS_AVG;
NOISE(NOISE > 650) = 0;

% Do IFFT to get noise in time domain
noise = ifft(NOISE);
noise(noise > 10) = 0;

% Plot
figure(6)
clf
hold on
stem(t, noise)
title('Task 3: noise(t)', 'FontSize', 12)
xlabel('t [sec]','FontSize',20);
ylabel('noise','FontSize',20);
xlim([0 length(p1_U_clean)*Ts])
ylim([-10 10])

% Calculate variance of the noise
noise_var = 0;
noise_var = var(noise);

% Explain
fprintf('In order to estimate the noise I take the averaged FFT of all five\nwindows and I set all the frequencies that belong to the real signal\nequal to zero. Those are the frequencies that were found in Task 2.\n') 
fprintf('Afterwards I do the IFFT of that signal, which only has the\nfrequency components of the noise, and calculate the variance of the\nnoise in the time domain.\n')
fprintf('\n')
fprintf('Figure 5 shows the reconstructed signal with no noise vs. the noisy\nsignal from Window 1. One can see that both signals are very similar\nappart from a small phase shift.\n')
fprintf('\n')
fprintf('Figure 6 shows the noise in the time domain.\n')
fprintf('\n')
fprintf('----------------------------------------------------------------------\n')
fprintf('----------------------------------------------------------------------\n')


p1_var = [noise_var];

end

function [autocorr_per, lags] = R_per(u, L)

autocorr_per = zeros(4*L,1);
for tau = (-2*L+1):(2*L)
    for i=0:L-1
            autocorr_per(tau+L*2) = autocorr_per(tau+L*2) + 1/L*u(i+1)*u(mod(i-tau,L)+1);
    end
lags = (-2*L+1):(2*L);
end
end

