% Differential evolution algorithm code, traverse ub and lb in Excel file to optimize and store the results.
clc;
clear all;
close all;

%% Step 1:Load ub and lb data in Excel file.
excel_data = readmatrix('Optimization_Range_0.1.xlsx');
[num_rows, num_cols] = size(excel_data);                         
assert(num_cols == 14, 'Each row in the Excel file should contain 14 values, of which the first 7 are lb and the last 7 are ub.');

%% Step 2:Load the data in all_cal_data_results.mat 
load('all_cal_data_results_Y_0.1.mat', 'all_results_Y');
assert(num_rows == size(all_results_Y, 1), 'The number of rows in Excel file does not match the number of rows in all_results.');
load('all_cal_data_results_X1_0.1.mat', 'all_results_X1');
load('all_cal_data_results_X2_0.1.mat', 'all_results_X2');


%% Step 3:Get input parameters
MaxIter = input('Maximum number of iterations of each data experiment£º');   %The maximum number of iterations for each independent experiment can be set to 100 or 200.
MC_Repeat = input('Number of independent experiments for each group of data£º');  %The number of times to run independent experiments can be set to 1 if the timeliness is considered.
fobj_template = @(r, Y0,X1,X2) FLNMGM_1N(r, Y0,X1,X2); 

%% Step 4:Initialization result storage
final_results = struct('final_best_score', [], 'final_best_position', [], 'ub', [], 'lb', []);

%% Step 5: Parameter optimization
for row_idx = 1:num_rows
    % Extract ub and lb of current row.
    lb = excel_data(row_idx, 1:7);
    ub = excel_data(row_idx, 8:14);

    % Extracting corresponding experimental data
    Y0 = all_results_Y(row_idx, :);
    X1=all_results_X1(row_idx, :);
    X2=all_results_X2(row_idx, :);
        
    fobj = @(r) fobj_template(r, Y0,X1,X2);

    pop = 60;
    dim = length(lb);
    all_results_scores = zeros(MC_Repeat, 1);
    all_positions = zeros(MC_Repeat, dim);

    % Run independent experiments
    for mc_run = 1:MC_Repeat
        fprintf('Optimization in line %d, running independently for the %d time.', row_idx, mc_run);
        [Best_Pos, Best_Score, IterCurve] = DE(pop, dim, ub, lb, fobj, MaxIter);
        all_results_scores(mc_run) = Best_Score;
        all_positions(mc_run, :) = Best_Pos;
    end

   % Statistical analysis results
    mean_score = mean(all_results_scores);
    std_score = std(all_results_scores);
    best_run = find(all_results_scores == min(all_results_scores), 1);
    final_best_position = all_positions(best_run, :);
    final_best_score = all_results_scores(best_run);

    % Save results
    final_results(row_idx).final_best_score = final_best_score;
    final_results(row_idx).final_best_position = final_best_position;
    final_results(row_idx).ub = ub;
    final_results(row_idx).lb = lb;

    % Display current results
    fprintf('Line %d is optimized, and the global optimal objective function value is achieved.: %.4f', row_idx, final_best_score);
end

%% Step 6:Save all results to a file
save('final_optimization_results.mat', 'final_results');
fprintf('All optimization is completed, and the results have been saved to final_optimization_results.mat');


%% DE optimization algorithm
function [Best_Pos, Best_Score, IterCurve] = DE(pop, dim, ub, lb, fobj, MaxIter)
 % Initializing population position    
    Positions = initialization(pop, dim, ub, lb);
    Best_Score = inf;
    Best_Pos = zeros(1, dim);
    IterCurve = zeros(1, MaxIter);
    F = 0.8; % Scaling factor
    CR = 0.9; % Crossover probability

    for t = 1:MaxIter
        for i = 1:pop
            % Select three random individuals
            idxs = randperm(pop, 3);
            while any(idxs == i)
                idxs = randperm(pop, 3);
            end
            x1 = Positions(idxs(1), :);
            x2 = Positions(idxs(2), :);
            x3 = Positions(idxs(3), :);

            % Mutation operation
            mutant = x1 + F * (x2 - x3);

             % Boundary check
            mutant = BoundrayCheck(mutant, ub, lb);

             % interlace operation
            trial = Positions(i, :);
            for j = 1:dim
                if rand() < CR
                    trial(j) = mutant(j);
                end
            end

             % Boundary check
            trial = BoundrayCheck(trial, ub, lb);

             % Selection operation
            if fobj(trial) < fobj(Positions(i, :))
                Positions(i, :) = trial;
            end

            % Update global optimal solution
            fitness = fobj(Positions(i, :));
            if fitness < Best_Score
                Best_Score = fitness;
                Best_Pos = Positions(i, :);
            end
        end
        IterCurve(t) = Best_Score;
        fprintf('Iteration times: %d, the current optimal fitness value: %.4f', t, Best_Score);
    end
end

function Positions = initialization(pop, dim, ub, lb)
    Positions = rand(pop, dim) .* (ub - lb) + lb;
end

function Position = BoundrayCheck(Position, ub, lb)
    Position = max(Position, lb);
    Position = min(Position, ub);
end
