function [p2_umin, p2_umax, p2_M, p2_Np, p2_u_etfe, p2_y_etfe, p2_omega, p2_estimate_no_window, p2_estimate_windowed, p2_gamma_best] = HS2019_SysID_midterm_p2_13921002()

    %% Solution for Problem 2
    
    %% General instructions for solution
    % Change the filename of this function, both in the function definition
    % above and in the filename in the folder

    % Modify your code in the next sections, and return the variables
    % requested.

    % If you skip one part of the problem, return the empty vectors as already
    % provided in the code

    % use the plant-function like this (where p2_u is any input you wish 
    % to apply and p2_y is the corresponding output of the plant):    
    % p2_y = HS2019_SysID_midterm_p2_system_sim(LegiNumber,p2_u);

    % use the validation-function like this (where p2_M is the period length 
    % of p2_u_etfe you found in Task 2):
    % p2_u_cv = HS2019_SysID_midterm_p2_validation(p2_M);

    % Extract Legi from Filename
    name = mfilename;
    LegiNumber = str2double(name(end-7:end));
    
    %% Task 1: Estimate plant input saturation 
    fprintf('\n')
    fprintf('----------------------------------------------------------------------\n')
    fprintf('----------------------------------------------------------------------\n')
    fprintf('\n')
    fprintf('Problem 2')
    fprintf('\n')
    fprintf('\n')
    fprintf('----------------------------------------------------------------------\n')
    fprintf('\n')
    fprintf('Task 1: Estimation of input saturation limits\n')
    fprintf('\n')
    
    Ts = 1;
    limits = -5:0.1:5;
    l_sig = 2000;
    
    p2_y = zeros(length(limits), l_sig);
    y_avg_vec = zeros(1,length(limits));
    for i = 1:101
        p2_u = limits(i) * ones(length(limits), l_sig);
        y_task1 = HS2019_SysID_midterm_p2_system_sim(LegiNumber,p2_u(i,:));
        y_avg_vec(1,i) = mean(y_task1(:,1));
    end
    
    % Plot
    figure(1)
    clf
    hold on
    title('Task 1: Estimation of input saturation limits', 'FontSize', 12)
    plot(limits, y_avg_vec)
    xlim([-5 5])
    ylim([-1.4 1.4])
    xlabel('amplitude u(t)', 'FontSize', 20)
    ylabel('average y(t)', 'FontSize', 20)
    
    umin = -1.1;
    umax = 1.2;
    
    % Explain
    fprintf('I used a step input signal of length 2000 samples with all\ndifferent possible magnitudes.\n')
    fprintf('Afterwards I plotted the average response and clearly saw what the\nsaturation limits are from the graph, in this case u_min = -1.1\nand u_max = 1.2. See Figure 1.\n')
    
    p2_umin = [-1.1]; %rounded to first decimal
    p2_umax = [1.2]; %rounded to first decimal
    
    
    %% Task 2: Design correct input
    fprintf('\n')
    fprintf('----------------------------------------------------------------------\n')
    fprintf('\n')
    fprintf('Task 2: Design of PRBS signal\n')
    fprintf('\n')
    
    M = 511; %from: 2*pi/L = pi/200, since for PRBS signals M=2^X-1
    
    u = 1.1 * idinput([M 1 7],'prbs');
    u_per = u(1:M);
    U_PER  = fft(u_per);
    U_PER_ABS = abs(U_PER);
    
    Np_theory = 0.25/(0.04/511*U_PER_ABS(2,1).^2);
    Np_theory = ceil(Np_theory);
    Np = Np_theory+1;
    
    y = HS2019_SysID_midterm_p2_system_sim(LegiNumber,u);
    
    % Explain
    fprintf('Constraint1: ')
    fprintf('Since the frequency resolution has to be larger or equal\npi/200 and we want our signal to be as short as possible, we follow\nthe formula in slide 2.18, so M has to be at least 400 samples.\n')
    fprintf('Since I have to design a PRBS signal and PRBS signals have a length\nof M=2^X-1, so the smallest number that fullfills both constraints\nis M=511.\n')
    fprintf('\n')
    
    fprintf('Constraint2: ')
    fprintf('I used the formulas in slides 4.4 and 4.3 to estimate\nthe minimal number of periods needed to achieve a variance of an\nunsmoothed transfer function estimate smaller than 0.04.\n')
    fprintf('From the formulas I get the minimal integer number of periods has to\nbe 6 and since I want to not regard the first period of the output to\nget rid of most of the transient, I take Np = 7.\n')
    fprintf('\n')
    
    fprintf('Constraint3: ')
    fprintf('From slide 3.15 we learned the ETFE alway gives an\nunbiased estimate for periodic input signals.\n')
    fprintf('\n')
    
    fprintf('Regarding all the constraints mentioned before, the input signal\nused is a PRBS signal with amplitude ±1.1 (since taking more would\ngo over the saturation limits), period length M = 511 and Np = 7.\n')
    
    p2_M       = [M]; % integer scalar, length of one period in p2_u_etfe
    p2_Np      = [Np_theory]; % integer scalar, number of periods in p2_u_etfe
    p2_u_etfe  = [u(M+1:length(u))]; % vector, your designed input signal used for the ETFE
    p2_y_etfe  = [y(M+1:length(y))]; % vector, the output from the plant, used for the ETFE


    %% Task 3: Compute unsmoothed ETFE
    fprintf('\n')
    fprintf('----------------------------------------------------------------------\n')
    fprintf('\n')
    fprintf('Task 3: Unsmoothed ETFE\n')
    fprintf('\n')
    
    t = 0:(length(u)-1);
    
    % Plot
    figure(2)
    clf
    hold on
    title('Task 3: Input and Output', 'FontSize', 12)
    plot(t, u)
    plot(t, y)
    legend('u(t)','y(t)')
    xlim([0 length(t)])
    ylim([-3 3])
    xlabel('t', 'FontSize', 20)
    ylabel('signal', 'FontSize', 20)
    
    % Throw away first period of y to get rid of transient and average to
    % get the noise. Then use one period to do the FFT.
    y_per_sum = zeros(M,1);
    for i = 2:Np
        for n = 1:M
            y_per_sum(n,1) = y_per_sum(n,1) + y((i-1)*M+n,1);
        end  
    end
    y_per_avg = y_per_sum / 5;
    
    Y_PER_AVG = fft(y_per_avg);
    
    % Get one period of u and do the FFT
    G_est = Y_PER_AVG ./ U_PER;
    
    omega = (2*pi/(Ts*M))*transpose(0:M-1);
    idx = find(omega > 0 & omega < pi);
    
    phase = zeros(511,1);
    phase = atan2(imag(G_est),real(G_est));
    
    % Plot G_est unsmoothed
    figure(3)
    clf
    hold on
    subplot(2,1,1);
    loglog(omega(idx),abs(G_est(idx)))
    title('Task 3: Bodeplot of ETFE G_e_s_t (unsmoothed)', 'Fontsize', 12)
    xlabel('\omega', 'FontSize', 20)
    ylabel('|G_u_n_s_m_o_o_t_h_e_d(e^{j\omega})|', 'FontSize', 20)
    xlim([0 pi])
    ylim([0 10])
    subplot(2,1,2);
    semilogx(omega(idx),phase(idx))
    xlabel('\omega', 'FontSize', 20)
    ylabel('\phi(e^{j\omega})', 'FontSize', 20)
    xlim([0 pi])
    ylim([-5 6])
    
    % Explain
    fprintf('I insert my 7-period-long PRBS signal in the\n"HS2019_SysID_midterm_p2_system_sim"-function to get the output of\nthe system. ')
    fprintf('I drop the first period of the output to get rid of the\ntransient. Afterwards, I average the 6 other periods of the output\nin order to get rid of noise.\n')
    fprintf('Then, I calculate the FFT of the averaged output, which is one period\nlong, and the FFT of one period of the input.\n')
    fprintf('Finally I divide the FFT of the averaged outoput by the FFT of one\nperiod of the input.\n')
    
    p2_omega               = [omega(idx)]; % vector, equally spaced frequecies in (0,pi)(not inclusive) 
    p2_estimate_no_window  = [G_est(idx)]; % vector, ETFE estimate unsmoothed
    
    
    %% Task 4: Compute smoothed ETFE with best Hann window width
    fprintf('\n')
    fprintf('----------------------------------------------------------------------\n')
    fprintf('\n')
    fprintf('Task 4: Smoothed ETFE with optimal γ for Hann window\n')
    fprintf('\n')
    
    % Create G_s with different gammas
    gamma = 50:50:500;
    
    % calculate M point FFTs
   
    u_cv = HS2019_SysID_midterm_p2_validation(p2_M);
    u_cv_per = u_cv(1:M);
    U_CV_PER = fft(u_cv_per);
    
    G_est = Y_PER_AVG ./ U_PER;% ETFE estimate
    Gs = zeros(length(gamma),length(G_est));% smoothed estimate
    
    y_cv = HS2019_SysID_midterm_p2_system_sim(LegiNumber,u_cv);
    
    y_cv_avg = zeros(M,1);
    y_cv_sum = zeros(M,1);
    for i = 1:5
        for n = 1:M
            y_cv_sum(n,1) = y_cv_sum(n,1) + y_cv(M*i+n,1);
        end
    end
    y_cv_avg = y_cv_sum / 5;
    
    Y_CV_PER = fft(y_cv_avg);
    
    % Create smoothed estimates
    for gam = 50:50:500
        [omega,Wg] = WfHann(gam,M);% window (centered)
        zidx = find(omega==0);% shift to start at zero
        omega = [omega(zidx:M);omega(1:zidx-1)];% frequency grid
        Wg = [transpose(Wg(zidx:M));transpose(Wg(1:zidx-1))];
        a = U_PER.*conj(U_PER);% variance weighting
    
        for wn = 1:M
            Wnorm = 0;% reset normalisation
                for xi = 1:M
                    widx = mod(xi-wn,M)+1;% wrap window index
                    Gs(gam/50,wn) = Gs(gam/50,wn) + Wg(widx) * G_est(xi) * a(xi);
                    Wnorm = Wnorm + Wg(widx) * a(xi);
                end
            Gs(gam/50,wn) = Gs(gam/50,wn)/Wnorm;% weight normalisation
        end
    end
    
    %Plot all smoothed estimates
    figure(4)
    clf
    hold on
    set(gca, 'XScale', 'log', 'YScale', 'log');
    for gam = 50:50:500
        loglog(omega(idx),abs(Gs(gam/50,idx)), 'DisplayName', ['Gamma:' num2str(gam)]);
    end
    legend
    title('Task 4: Smoothed ETFE G_e_s_t with different Hann window sizes', 'Fontsize', 12)
    xlabel('\omega', 'FontSize', 20)
    ylabel('|G_e_s_t|', 'FontSize', 20)
    xlim([0 pi])
    ylim([0 10])
    
    error = zeros(10,1);
    
    for gam = 1:10
            error(gam,1) = norm((Y_CV_PER - transpose(Gs(gam,:)) .* U_CV_PER));
    end
        
    % Plot errors
    figure(5)
    clf
    hold on
    title('Task 4: Errors for the different window sizes γ', 'FontSize', 12)
    plot(50:50:500, abs(error))
    xlim([50 500])
    ylim([118 136])
    xlabel('γ', 'FontSize', 20)
    ylabel('error', 'FontSize', 20)
    
    % Explain
    fprintf('The error changes with γ following the bias-variance trade-off\ndescribed in slide 4.12. ')
    fprintf('In the first part of the curve, the MSE\n(here SE) goes down because of the decay of the variance (1/Np)\nuntil it reaches its minimum. After that point, the quadratic\ninfluence of the bias (grlows linearly) will make the MSE rise.\n')
    fprintf('\n')
    fprintf('If I was taking more periods for my PRBS input signal, the variance\nwould be smaller due to smaller noise, so the optimal Hann window\nwould have a smaller optimal γ-value.\n')
    fprintf('\n')
    fprintf('Figure 4 shows the smoothed ETFE for different γ.\n')
    fprintf('\n')
    fprintf('Figure 5 shows the error function when smoothing G with different γ.\n')
    fprintf('\n')
    fprintf('----------------------------------------------------------------------\n')
    fprintf('----------------------------------------------------------------------\n')
    
    p2_estimate_windowed = [Gs(4,:)]; %vector, Hann window smoothed estimate using window width p2_gamma_best
    p2_gamma_best        = [200]; %scalar integer, best Hann window width


end







