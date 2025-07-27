%% Generation of experimental data with added noise

%% Step 1:Experimental parameter setting
num_experiments = 200;  % Set the number of experiments
f=[0.5 10];  %Nakagami function parameters
Y0(1,1)=1;   %Initial value of dependent variable
Row2=2;      %Number of independent variables
S=[1,4];
r=[0.5,1.5,0.8,1.1]; %Lag parameter and power index value
miu_1=0.9;  
miu_2N=[0.5,0.7];
miu_3=0.2;
miu_4=0.3;
len=30;   %Dependent variable data length
a=1;    %Noise standard deviation
all_results = zeros(num_experiments, len); % Store the results of all experiments.

%% Step 2: Weight function matrix
M=[];
N=[];
const = 2 * f(1)^f(1) / (gamma(f(1)) * f(2)^f(1));  % Unified coefficient constant
for i = 1:Row2
    ri = abs(r(i));
    ti = floor(ri);
    t(i) = ti;

    j = 0:ti;
    % Calculate the weight corresponding to item j.
    tmp1 = const * (j.^(2*f(1)-1)) .* exp(-f(1) * j.^2 / f(2));
    tmp2 = sum(tmp1);
    M(i,1:ti+1) = tmp1;  

    if ri == ti  % The lag parameter is an integer
        N(i) = tmp2;
    else  % Lag parameter is non-integer
        tmp3 = const * (ri.^(2*f(1)-1)) * exp(-f(1) * ri.^2 / f(2));
        M(i,ti+2) = tmp3;
        N(i) = tmp2 + tmp3;
    end
end


Beta=[];
for i = 1:Row2
    t(i) = floor(abs(r(i)));
    
    for j = 0:t(i)
        Beta(i, j+1) = M(i, j+1) / N(i);
    end

    % Lag parameter is non-integer
    if r(i) ~= floor(r(i))
        Beta(i, t(i)+2) = M(i, t(i)+2) / N(i);
    end
end

%% Step3: Dependent variable and independent variable experimental data generation
for exp_idx = 1:num_experiments
    % Initializing random error term
    epsilonx1 = normrnd(0, sqrt(a), len+1, 1);
    epsilonx2 = normrnd(0, sqrt(a), 4*len+4, 1);
    epsilony = normrnd(0, sqrt(a), len, 1);
 
        x1_0_n = zeros(1, len+1); % Initialize the result vector
        x2_0_n = zeros(1, 4*len+4); % Initialize the result vector
    for i = 1:len+1
        x1_0_n(i) = 0.6 * exp(0.15 * i) + epsilonx1(i);
    end
    for i = 1:4*len+4
        x2_0_n(i) = 0.5 * i + 2 * sin(2 * pi * (i / 4)) + 0.8 * sin(2 * pi * (i / 8)) + epsilonx2(i);
    end
   
    maxLength =length(x2_0_n);

    % Fill a shorter vector with 0.
    X1_padded = [x1_0_n, zeros(1, maxLength - length(x1_0_n))];  
    X2_padded = [x2_0_n, zeros(1, maxLength - length(x2_0_n))];  

    % Splice into a matrix
    Xs = [X1_padded; X2_padded];
    Xsr_AGO=cumsum(Xs, 2);
      
    
    % Initialize Simu_Ys_r
    Simu_Ys_r = zeros(1, len);
    Simu_Ys_r(1) = Y0(1, 1);

    % Generate Simu_Ys_r data.
    for k = 2:len
        tmp1 = power(miu_1, k-1) * Y0(1, 1);
        tmp2 = 0;
        for w = 2:k
            tmp2 = tmp2 + power(miu_1, k-w) * (miu_3 * w + miu_4);
        end
        tmp6 = 0;
        for v = 2:k
            tmp5 = 0;
            for i = 1:Row2
                t(i) = floor(abs(r(i)));
                tmp3 = 0;
                if r(i) >= 0
                    for j = 0:t(i)
                        tmp3 = tmp3 + Beta(i, j+1) * Xsr_AGO(i, v*S(i)-j);
                    end
                    if abs(r(i)) == t(i)
                        tmp4 = tmp3;
                    else
                        tmp4_1 = tmp3 + Beta(i, t(i)+2) * ((Xsr_AGO(i, k*S(i)-floor(r(i))) - Xsr_AGO(i, k*S(i)-(floor(r(i))+1))) * power(floor(r(i))+1-r(i), 1) + Xsr_AGO(i, k*S(i)-(floor(r(i))+1)));
                        tmp4 = power(tmp4_1, r(Row2+i));
                    end
                    tmp5 = tmp5 + tmp4 * miu_2N(i);
                else
                    for j = 0:t(i)
                        tmp3 = tmp3 + Beta(i, j+1) * Xsr_AGO(i, v*S(i)+j);
                    end
                    if abs(r(i)) == t(i)
                        tmp4 = tmp3;
                    else
                        tmp4_1 = tmp3 + Beta(i, t(i)+2) * ((Xsr_AGO(i, k*S(i)-floor(r(i))) - Xsr_AGO(i, k*S(i)-(floor(r(i))+1))) * power(floor(r(i))+1-r(i), 1) + Xsr_AGO(i, k*S(i)-(floor(r(i))+1)));
                        tmp4 = power(tmp4_1, r(Row2+i));
                    end
                    tmp5 = tmp5 + tmp4 * miu_2N(i);
                end
            end
            tmp6 = tmp6 + power(miu_1, k-v) * tmp5;
        end
        tmp = tmp1 + tmp2 + tmp6 + epsilony(k);
        Simu_Ys_r(k) = tmp * Y0(1, 1);
    end

    All_cal_data = [Simu_Ys_r(1), diff(Simu_Ys_r)];
    all_results_Y(exp_idx, :) = All_cal_data;
    all_results_X1(exp_idx, :) = X1_padded;
    all_results_X2(exp_idx, :) = X2_padded;

end

%% Step 4: Save the results to a file
save('all_cal_data_results_Y.mat', 'all_results_Y');
save('all_cal_data_results_X1.mat', 'all_results_X1');
save('all_cal_data_results_X2.mat', 'all_results_X2');

%% Step 5: Show partial results
disp('Ç°5×é All_cal_data:');
disp(all_results_Y(1:5, :));
