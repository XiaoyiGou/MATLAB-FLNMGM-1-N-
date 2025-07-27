%% Traverse all_results array in all_cal_data_results.mat file line by line and perform phase spectrum analysis and coherence calculation.
clc, clear all;

%% Step1: Load data
load('all_cal_data_results_Y_1.mat', 'all_results_Y');
% load('all_cal_data_results_X_1.mat', 'all_results_X1');
load('all_cal_data_results_X2_1.mat', 'all_results_X2');
[num_rows, ~] = size(all_results_Y);

results = struct('max_lag_time', [], 'min_lag_time', []);  %Store results

for row_idx = 1:num_rows
%% Step 2: Set dependent and independent variables
    len=30;                         %Dependent variable data length
    Y = all_results_Y(row_idx, :);  %Dependent variable 

% %  Independent variable X1
%    X= all_results_X1(row_idx, 1:len);  
%    F=1;                                %Frequency multiple difference between dependent variable and independent variable

% % Independent variable X2
    F=4;                                %Frequency multiple difference between dependent variable and independent variable
    X2= all_results_X2(row_idx, 1:F*len);  
    % Make sure that the number of columns is a multiple of the frequency difference.
    num_cols = length(X2);
    if mod(num_cols, F) ~= 0
        error('The number of columns must be a multiple of F.');
    end
    X2_reshaped = reshape(X2, F, []);
    X = sum(X2_reshaped, 1); 
    X = X(1:len);     % Ensure that it is consistent with the dependent variable data length.
    
    
    % Ensure data consistency
    N = length(X);
    if length(Y) ~= N
        error('The data length of variable y and variable x must be equal£¡');
    end

    %% Step 3: Normalization processing
    X_norm = (X - mean(X)) / std(X);
    Y_norm = (Y - mean(Y)) / std(Y);

    %% Step 4: Computational Fourier transform
    X_fft = fft(X_norm);
    Y_fft = fft(Y_norm);

    %% Step 5:Set frequency axis
    fs = 1;                      % The data sampling interval is 1 year.
    frequencies = (0:N-1)*(fs/N);

    %% Step 6: Calculation of phase difference, lag time and coherence
    phase_diff = angle(Y_fft) - angle(X_fft);

    lag_time = phase_diff ./ (2 * pi * frequencies);
    lag_time(frequencies == 0) = NaN; % Avoid zero denominator

    [Cxy, f] = mscohere(Y_norm, X_norm, [], [], N, fs);

    %% Step 7:Extract the results in the range of 0-0.3Hz.
    lag_time_range = lag_time(frequencies <= 0.3);
    frequencies_range = frequencies(frequencies <= 0.3);
    Cxy_range = Cxy(f <= 0.3);
    f_range = f(f <= 0.3);

    max_lag_time = max(lag_time_range)*F;
    min_lag_time = min(lag_time_range)*F;
    
    if max_lag_time > 2*F - 1
        max_lag_time = 2*F - 1;
    end

    if min_lag_time < -(2*F - 1)
        min_lag_time = -(2*F - 1);
    end
    
    %% Step 8:Save the results of the current row.
    results(row_idx).max_lag_time = max_lag_time;
    results(row_idx).min_lag_time = min_lag_time;

    %% Step 9:Display results
    fprintf('Row %d: Max Lag Time (0-0.3Hz): %f years\n', row_idx, max_lag_time);
    fprintf('Row %d: Min Lag Time (0-0.3Hz): %f years\n', row_idx, min_lag_time);
end

%% Step 10:Save the results to a file
% save('phase_coherence_results_X1.mat', 'results');
save('phase_coherence_results_X2.mat', 'results');

fprintf('The analysis of all lines is completed, and the results have been saved in phase_coherence_results.mat\n');
