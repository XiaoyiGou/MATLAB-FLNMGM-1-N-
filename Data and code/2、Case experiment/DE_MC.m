clc;
clear all;
close all;

%Jiangsu
pop = 30;
dim = 9;
ub =[-0.1250 -0.1396 5.2042 2 2 2 5 10 10];
lb =[-0.9695 -1.0416 -3.6028 0.0001 0.0001 0.0001 0.5 0];

%Guizhou
% pop = 30;
% dim = 9;
% ub =[-0.1502 -0.0812 23.0000 2 2 2 5 10 10];
% lb =[-0.8822 -0.7181  -7.3405 0.0001 0.0001 0.0001 0.0001 0.5 0];

MaxIter = input('Maximum number of iterations£º');   %The maximum number of iterations for each independent experiment can be set to 100 or 200.
MC_Repeat = input('Independent operation times£º');  %The number of times to run independent experiments can be set to 1 if the timeliness is considered.
fobj = @(r) FLNMGM_1N(r); 

% Start the time
tic;
all_results = zeros(MC_Repeat, 1);
all_positions = zeros(MC_Repeat, dim);

for mc_run = 1:MC_Repeat
    fprintf('Run independently for the %d time\n', mc_run);
    [Best_Pos, Best_Score, IterCurve] = DE(pop, dim, ub, lb, fobj, MaxIter);
    all_results(mc_run) = Best_Score;
    all_positions(mc_run, :) = Best_Pos;
end

elapsedTime = toc;

% Statistical analysis results
mean_score = mean(all_results);
std_score = std(all_results);
best_run = find(all_results == min(all_results), 1);
final_best_position = all_positions(best_run, :);
final_best_score = all_results(best_run);

fprintf('The code running time is %.4f seconds\n', elapsedTime);
fprintf('Average objective function value: %.4f\n', mean_score);
fprintf('Standard deviation of objective function value: %.4f\n', std_score);
fprintf('Global optimal objective function value: %.4f\n', final_best_score);
fprintf('Corresponding optimal solution: %s\n', mat2str(final_best_position));


figure(1);
plot(IterCurve, 'r-', 'linewidth', 1);
grid on;
title('Iterative curve of differential evolution algorithm');
xlabel('Iterations');
ylabel('Fitness value');

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
        fprintf('Iteration times: %d, the current optimal fitness value: %.4f\n', t, Best_Score);
    end
end

function Positions = initialization(pop, dim, ub, lb)
    Positions = rand(pop, dim) .* (ub - lb) + lb;
end

function Position = BoundrayCheck(Position, ub, lb)
    Position = max(Position, lb);
    Position = min(Position, ub);
end
