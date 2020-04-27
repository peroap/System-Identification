function [p1_theta_est1,p1_Phi,p1_theta_est2,p1_y_pred] = HS2019_SysID_final_p1_13921002()
%% Solution for Problem 1
%% Output format specification
% p1_R must be a 1xT vector
% p1_omega must be a 1xM vector
% p1_a must be a 1xM vector
% p1_var must be a scalar
%% Generate data

% Extract Legi from Filename
name=mfilename;
LegiNumber= str2num(name(end-7:end));

[p1_u, p1_y, p1_theta_hat, p1_u_past, p1_y_past,p1_pred_err] = HS2019_SysID_final_p1_GenerateData(LegiNumber);

%% General instructions for solution

% Change the filename of this function, both in the function definition
% above and in the filename in the folder

% Use the variable p1_U to solve the problem. 

% Modify your code in the next sections, and return the variables
% requested.

% If you skip one part of the problem, return the empty vectors as already
% provided in the code

%% Task 1: Obtain initial estimate

fprintf("------------------------------------------------------------------\n")
fprintf("------------------------------------------------------------------\n")
fprintf("------------------------------------------------------------------\n")
fprintf("\n")
fprintf("*** Problem 1 ***\n")
fprintf("\n")

fprintf("------------------------------------------------------------------\n")
fprintf("------------------------------------------------------------------\n")
fprintf("\n")
fprintf("Task 1:\n")
fprintf("\n")

fprintf("Since I can assume c_1 and c_2 are zero, the model I get is an ARX\nmodel.\n")
fprintf("Regarding the fact that there is correlated noise in the system,\n")
fprintf("the best approach is to use the Markov Estimator (also called BLUE)\n")
fprintf("as in slide 9.32.\n")
fprintf("Since A is monic, I transform the equation to get the Phi-matrix.\n")
fprintf("The R-matrix is built following the definition of the noise\ncorrelation.\n")

fprintf("\n")
fprintf("Equations for the computation of theta_est1:\n")
fprintf("theta_est1 = transpose(Zeta1) * y;\n")
fprintf("where:\n")
fprintf("Zeta1 = inv(R) * Phi * inv(transpose(Phi) * inv(R) * Phi), as in 9.32\n")

fprintf("\n")
fprintf("theta_est1 is biased, motivated from the assumption that the system\n")
fprintf("is linear and ARX, when in reality it's not. A similar example was\n")
fprintf("introduced during lecture 9 in sl 9.36.\n")
fprintf("The covariance of theta_est is smaller or equal than the covariance\n")
fprintf("of any other estimator, as stated in 9.32 (Markov Estimator).\n")
fprintf("\n")

% Markov BLUE estimator

N = length(p1_u);

R = zeros(N);
for i=1:N
    for k=1:N
        if abs(i-k)<4
            R(i,k) = 0.2;
        end
        if i == k
            R(i,k) = 0.9;
        end
    end
end

theta_ls = zeros(6,1);

Phiy_p1 = zeros(N,6);
Phiy_p1(1,:) = [0,0,0,0,0,0];
Phiy_p1(2,:) = [-p1_y(1),0,0,p1_u(1),0,0];
Phiy_p1(3,:) = [-p1_y(2),-p1_y(1),0,p1_u(2),p1_u(1),0];
for i = 4:N
    Phiy_p1(i,:) = [-p1_y(i-1),-p1_y(i-2),-p1_y(i-3),p1_u(i-1),p1_u(i-2),p1_u(i-3)];
end

theta_ls(:,1) = Phiy_p1\p1_y;


zeta_p1 = inv(R) * Phiy_p1 * inv(transpose(Phiy_p1) * inv(R) * Phiy_p1);

theta_blue_1 = zeros(6,1);
theta_blue_1 = transpose(zeta_p1) * p1_y;



p1_theta_est1 = [theta_blue_1];
p1_Phi = [Phiy_p1];

%% Task 2: Improve your estimate

fprintf("------------------------------------------------------------------\n")
fprintf("------------------------------------------------------------------\n")
fprintf("\n")
fprintf("Task 2:\n")
fprintf("\n")

fprintf("Since the coefficients c_1 and c_2 are given, the ARMAX system can\n")
fprintf("be transformed into a ARX system. In order to have a 'truly linear'\n")
fprintf("system, I have to transform input u and output y to fit the new\n")
fprintf("ARX model.\n")
fprintf("\n")

fprintf("Equations for the computation of theta_est2:\n")
fprintf("The equation for the new system is:\n")
fprintf("inv(C) * y = B/A * inv(C) * u + 1 A(z) A(z)e(k)\n")
fprintf("\n")
fprintf("One can see that if\n")
fprintf("y_new = inv(C) * y\n")
fprintf("and\n")
fprintf("u_new = inv(C) * u\n")
fprintf("are defined, the system has an ARX structure (as in ExerciseSet8).\n")
fprintf("Zeta2 = inv(R) * Phi_new * inv(transpose(Phi_new) * inv(R) * Phi_new)\n")
fprintf("\n")
fprintf("The estimator is computed like this:\n")
fprintf("theta_est2 = transpose(Zeta2) * y_new;\n")
fprintf("\n")
fprintf("The estimator is in this case unbiased, since the assumption of a\n")
fprintf("linear system holds, because of the transformation of u and y.\n")


%Define u_new and y_new as in ExSet8 solution

c_1 = 1;
c_2 = pi/4;

C = [1 c_1 c_2];
C_inv = tf(1,C,1,'Variable','z^-1');

u_new = lsim(C_inv, p1_u);
y_new = lsim(C_inv, p1_y);

%Define new Phi
Phiy_p2 = zeros(N,6);
Phiy_p2(1,:) = [0,0,0,0,0,0];
Phiy_p2(2,:) = [-y_new(1),0,0,u_new(1),0,0];
Phiy_p2(3,:) = [-y_new(2),-y_new(1),0,u_new(2),u_new(1),0];
for i = 4:N
    Phiy_p2(i,:) = [-y_new(i-1),-y_new(i-2),-y_new(i-3),u_new(i-1),u_new(i-2),u_new(i-3)];
end

%Find Markov Estimator with new Phi

zeta_p2 = inv(R) * Phiy_p2 * inv(transpose(Phiy_p2) * inv(R) * Phiy_p2);

theta_blue_2 = zeros(6,1);
theta_blue_2 = transpose(zeta_p2) * y_new;

theta_2_ls = Phiy_p2 \ y_new;

p1_theta_est2 = [theta_blue_2];


%% Task 3: Compute prediction

fprintf("------------------------------------------------------------------\n")
fprintf("------------------------------------------------------------------\n")
fprintf("\n")
fprintf("Task 3:\n")
fprintf("\n")

fprintf("I use the formula to perform the one-step ahead prediction (10.17).\n")
fprintf("Since the prediction error in this case is dependent of theta,\n")
fprintf("the predictor was calculated componentwise as in 10.18:\n")
fprintf("y_est_3 = b_1*p1_u_past(1) + b_2*p1_u_past(2) + b_3*p1_u_past(3) + ...\n")
fprintf("(-a_1)*p1_y_past(1) + (-a_2)*p1_y_past(2) + (-a_3)*p1_y_past(3) + ...\n")
fprintf("+ c_1*p1_pred_err(1) + c_2*p1_pred_err(2)\n")

%Create A_est and B_est with estimated parameters from part 2

a_1 = p1_theta_hat(1,1);
a_2 = p1_theta_hat(2,1);
a_3 = p1_theta_hat(3,1);
b_1 = p1_theta_hat(4,1);
b_2 = p1_theta_hat(5,1);
b_3 = p1_theta_hat(6,1);

z = tf('z');
A_est = 1 + a_1*z^(-1) + a_2*z^(-2) + a_3*z^(-3);
B_est = b_1*z^(-1) + b_2*z^(-2) + b_3*z^(-3);

y_est_3 = b_1*p1_u_past(1) + b_2*p1_u_past(2) + b_3*p1_u_past(3) + (-a_1)*p1_y_past(1) ...
+ (-a_2)*p1_y_past(2) + (-a_3)*p1_y_past(3) + c_1*p1_pred_err(1) + c_2*p1_pred_err(2);



p1_y_pred = [y_est_3];
end
