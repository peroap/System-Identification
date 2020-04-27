function [p3_b_ML,p3_b_MAP,p3_cv_error,p3_prior_best] = HS2019_SysID_final_p3_13921002()

    %% Solution for Problem 3

    %% General instructions for solution
    % Change the filename of this function, both in the function definition
    % above and in the filename in the folder

    % Modify your code in the next sections, and return the variables
    % requested.

    % If you skip one part of the problem, return the null variables as already
    % provided in the code

    % Extract Legi from Filename
    name = mfilename;
    LegiNumber = str2double(name(end-7:end));
    
    % Obtain experiment data
    [p3_u,p3_y,p3_u_cv,p3_y_cv] = HS2019_SysID_final_p3_GenerateData(LegiNumber);

    %% Task 1: Maximum likelihood estimate
    fprintf("------------------------------------------------------------------\n")
    fprintf("------------------------------------------------------------------\n")
    fprintf("------------------------------------------------------------------\n")
    fprintf("\n")
    fprintf("*** Problem 3 ***\n")
    fprintf("\n")
    
    fprintf("------------------------------------------------------------------\n")
    fprintf("------------------------------------------------------------------\n")
    fprintf("\n")
    fprintf("Task a:\n")
    fprintf("\n")
    
    fprintf("The likelihood function describes the probability of b_ML\n")
    fprintf("taking specific values. The MLE takes the most probable values\n")
    fprintf("of the likelihood function as the estimated values.\n")    
    fprintf("\n")
    fprintf("Equation for computing b_ML:\n")
    fprintf("b_ML = inv(transpose(Phi)*Phi)*transpose(Phi)*p3_y\n")
    fprintf("\n")
    fprintf("In order to compute the estimate, use the fact that\n")
    fprintf("y - sum(b(i)*u(k.i) ~ Normal(0,0.5^2)\n")
    fprintf("Because of that, I can calculate the Phi-matrix with 8 input\n")
    fprintf("time-shifts (since 8 b-parameters) and use the least squares\n")
    fprintf("approach to get b_ML (equation from above).\n")    
    fprintf("\n")


    N = length(p3_u);
    
    Phi = zeros(N,8);
    Phi(1,:) = [0,0,0,0,0,0,0,0];
    Phi(2,:) = [p3_u(1),0,0,0,0,0,0,0];
    Phi(3,:) = [p3_u(2),p3_u(1),0,0,0,0,0,0];
    Phi(4,:) = [p3_u(3),p3_u(2),p3_u(1),0,0,0,0,0];
    Phi(5,:) = [p3_u(4),p3_u(3),p3_u(2),p3_u(1),0,0,0,0];
    Phi(6,:) = [p3_u(5),p3_u(4),p3_u(3),p3_u(2),p3_u(1),0,0,0];
    Phi(7,:) = [p3_u(6),p3_u(5),p3_u(4),p3_u(3),p3_u(2),p3_u(1),0,0];
    Phi(8,:) = [p3_u(7),p3_u(6),p3_u(5),p3_u(4),p3_u(3),p3_u(2),p3_u(1),0];
    for i = 9:N
        Phi(i,:) = [p3_u(i-1),p3_u(i-2),p3_u(i-3),p3_u(i-4),p3_u(i-5),p3_u(i-6),p3_u(i-7),p3_u(i-8)];
    end
        
    theta_mle = zeros(8,1);
    theta_mle = inv(transpose(Phi)*Phi) * transpose(Phi)*p3_y;
   
    p3_b_ML = theta_mle;   % vector of dimension 8x1


    %% Task 2: Maximum a posteriori estimates
    fprintf("------------------------------------------------------------------\n")
    fprintf("------------------------------------------------------------------\n")
    fprintf("\n")
    fprintf("Task b:\n")
    fprintf("\n")
    
    
    %Following definition from ExSet10 Problem2
    fprintf("Optimization problem: Assuming a prior knowledge of the\n")
    fprintf("estimator (S), what will the best estimate be?\n")
    fprintf("As for b_ML also exploit the fact that\n")
    fprintf("y - sum(b(i)*u(k.i) ~ Normal(0,0.5^2);\n")
    fprintf("but including Bayes' rule. One gets:\n")
    fprintf("b_MAP = argmax_b(p(b|y)) = argmin_b(V(b));\n")
    fprintf("where\n")
    fprintf("V(b) := likelihood function of b.\n")
    fprintf("\n")
    fprintf("Equation for computing b_MAP:\n")
    fprintf("b_MAP=inv(transpose(Phi)*inv(sigma)*Phi+inv(S))*transpose(Phi)*inv(sigma)*y;\n")
    fprintf("where:\n")
    fprintf("sigma := covariance matrix of the data distribution\n")
    fprintf("S := prior knowledge about our estimators\n")
    fprintf("\n")
    
    fprintf("The difference between b_ML and b_MAP is that when calculating\n")
    fprintf("b_ML one only takes into account the likelihood function and\n")
    fprintf("takes the values with the biggest probability, and for computing\n")
    fprintf("b_MAP the prior knowledge is included via Bayes' rule. MAP 'biases'\n")
    fprintf("the estimate that one gets from only looking at the likelihood\n")
    fprintf("function towards the prior knowledge distribution.\n")
    fprintf("\n")
    
    fprintf("The different choices of S assume different distributions of the\n")
    fprintf("estimators in b_MAP. The estimators that has the lowest\n")
    fprintf("validation MSE is S_2, with a diagonal matrix (iid estimators,\n")
    fprintf("since zero-entries for i!=j).\n")

    sigma = 0.5^2 * eye(N);
        
    S_1 = eye(8);
    %Define S_2
    S_2 = zeros(8);
    for i = 1:8
        S_2(i,i) = 0.8^i;
    end
    %Define S_3
    S_3 = zeros(8);
    for i = 1:8
        S_3(i,i) = 0.5^i;
    end
    %Define S_4
    S_4 = zeros(8);
    for i = 1:8
        for j = 1:8
            if i >= j
                S_4(i,j) = 0.8^i;
            end
            if i < j
                S_4(i,j) = 0.8^j;
            end
        end
    end
    %Define S_5
    S_5 = zeros(8);
    for i = 1:8
        for j = 1:8
            if i >= j
                S_5(i,j) = 0.5^i;
            end
            if i < j
                S_5(i,j) = 0.5^j;
            end
        end
    end
    
    theta_map1 = inv(transpose(Phi) * inv(sigma) * Phi + inv(S_1))*transpose(Phi)*inv(sigma)*p3_y;
    theta_map2 = inv(transpose(Phi) * inv(sigma) * Phi + inv(S_2))*transpose(Phi)*inv(sigma)*p3_y;
    theta_map3 = inv(transpose(Phi) * inv(sigma) * Phi + inv(S_3))*transpose(Phi)*inv(sigma)*p3_y;
    theta_map4 = inv(transpose(Phi) * inv(sigma) * Phi + inv(S_4))*transpose(Phi)*inv(sigma)*p3_y;
    theta_map5 = inv(transpose(Phi) * inv(sigma) * Phi + inv(S_5))*transpose(Phi)*inv(sigma)*p3_y;
    
    b_map = zeros(8,5);
    b_map = [theta_map1 theta_map2 theta_map3 theta_map4 theta_map5];
    
    %Define y_ucv
    B_map1 = tf([0 theta_map1(1) theta_map1(2) theta_map1(3) theta_map1(4) theta_map1(5) theta_map1(6) theta_map1(7) theta_map1(8)],1,1,'Variable', 'z^-1');
    B_map2 = tf([0 theta_map2(1) theta_map2(2) theta_map2(3) theta_map2(4) theta_map2(5) theta_map2(6) theta_map2(7) theta_map2(8)],1,1,'Variable', 'z^-1');
    B_map3 = tf([0 theta_map3(1) theta_map3(2) theta_map3(3) theta_map3(4) theta_map3(5) theta_map3(6) theta_map3(7) theta_map3(8)],1,1,'Variable', 'z^-1');
    B_map4 = tf([0 theta_map4(1) theta_map4(2) theta_map4(3) theta_map4(4) theta_map4(5) theta_map4(6) theta_map4(7) theta_map4(8)],1,1,'Variable', 'z^-1');
    B_map5 = tf([0 theta_map5(1) theta_map5(2) theta_map5(3) theta_map5(4) theta_map5(5) theta_map5(6) theta_map5(7) theta_map5(8)],1,1,'Variable', 'z^-1');

    y_map1 = lsim(B_map1, p3_u_cv);
    y_map2 = lsim(B_map2, p3_u_cv);
    y_map3 = lsim(B_map3, p3_u_cv);
    y_map4 = lsim(B_map4, p3_u_cv);
    y_map5 = lsim(B_map5, p3_u_cv);
    
    mse1 = sum(1/N*(p3_y_cv - y_map1).^2);
    mse2 = sum(1/N*(p3_y_cv - y_map2).^2);
    mse3 = sum(1/N*(p3_y_cv - y_map3).^2);
    mse4 = sum(1/N*(p3_y_cv - y_map4).^2);
    mse5 = sum(1/N*(p3_y_cv - y_map5).^2);
    
    cv_error = zeros(5,1);
    cv_error(1,1) = mse1;
    cv_error(2,1) = mse2;
    cv_error(3,1) = mse3;
    cv_error(4,1) = mse4;
    cv_error(5,1) = mse5;
    
    p3_b_MAP        = b_map;   % matrix of dimension 8x5
    p3_cv_error    	= cv_error;   % vector of dimension 5x1
    p3_prior_best   = 2;            % scalar integer in the set {1,2,3,4,5}


end
