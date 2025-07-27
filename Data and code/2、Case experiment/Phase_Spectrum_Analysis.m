%Phase spectrum analysis, coherence

clc;
clear all;
close all;
%disp[--------------------------------------------¡¾Jiangsu¡¿--------------------------------------------]
% % %EC
Y=[5012.54	5114.70	5458.95	5807.89	6128.27	6264.36	6373.71	7101.16	7399.55	7832.96	8487.00];
% % %Population
X=[8281.09	8315.11	8381.47	8423.5	8446.19	8469.09	8477.26	8505.4	8515	8526	8526];
F=1;    % Converted back to the statistical frequency corresponding to the independent variable, and the lag period unit: year.
% % %GDP
% X=[65088.32	70116.38	76086.17	85900.94	93207.55	98656.82	102807.68	117392	122089.3	128222.2	137008];
% F=4;  % Converted back to the statistical frequency corresponding to the independent variable, and the lag period unit: quarter.
% % %Import and export trade
% X=[346642990	338966597	340902249.6	400333030	438788960.1	435561207	446066419.6	522356820.3	546458782.8	527081267.4	562640243.3];
% F=12; % Converted back to the statistical frequency corresponding to the independent variable, and the lag period unit: month.

%disp[--------------------------------------------¡¾Guizhou¡¿--------------------------------------------]
% % %EC
% Y=[1173.74	1174.21	1241.78	1384.89	1482.12	1540.68	1586	1743	1743	1783	1902];
% % % %Population
% X=[3677	3708	3758	3803	3822	3848	3858	3852	3856	3865	3860];
% F=1;     % Converted back to the statistical frequency corresponding to the independent variable, and the lag period unit: year.
% % % %GDP
% X=[9251.01	10502.56	11734.43	13540.83	15353.21	16769.34	17860.41	19459	20010.4	20913.25	22667.12];
% F=4;   % Converted back to the statistical frequency corresponding to the independent variable, and the lag period unit: quarter.
% % %Import and export trade
% X=[6652054	7919232	4830310.2	5536931.56	5040637.97	4522623	5521600.27	6541054.28	8027256.77	8194239.54	8665812.36];
% F=12;  % Converted back to the statistical frequency corresponding to the independent variable, and the lag period unit: month.


% Ensure data consistency
N = length(X);
if length(Y) ~= N
    error('The data length of variable y and variable x must be equal£¡');
end

% Normalization (average removal and scaling to the same range)
X = (X - mean(X)) / std(X);
Y = (Y - mean(Y)) / std(Y);

% Computational Fourier transform
X_fft = fft(X);
Y_fft = fft(Y);

% Frequency axis
fs = 1; % The data sampling interval is 1 year.
frequencies = (0:N-1)*(fs/N);

%Phase difference calculation
phase_diff = angle(Y_fft)-angle(X_fft);

%Lag time calculation
lag_time = phase_diff ./ (2 * pi * frequencies)*F; %Restore the corresponding statistical frequency
lag_time(frequencies == 0) = NaN; % Avoid zero denominator

% Computational coherence
[Cxy, f] = mscohere(Y, X, [], [], N, fs);

% Display 0-0.3Hz range
lag_time_range = lag_time(frequencies <= 0.3);
frequencies_range = frequencies(frequencies <= 0.3);
Cxy_range = Cxy(f <= 0.3);
f_range = f(f <= 0.3);

%Calculate the maximum and minimum values of the time range after the gradual change
max_lag_time = max(lag_time_range);
min_lag_time = min(lag_time_range);

if max_lag_time > 2*F - 1
    max_lag_time = 2*F - 1;
end

if min_lag_time < -(2*F - 1)
    min_lag_time = -(2*F - 1);
end


% Display results
fprintf('Max Lag Time (0-0.3Hz): %.4f \n', max_lag_time);
fprintf('Min Lag Time (0-0.3Hz): %.4f \n', min_lag_time);

% Visualization
figure;
% Phase difference diagram
subplot(3, 1, 1);
plot(frequencies, phase_diff, 'LineWidth', 1.5);
set(gca,'FontSize',20, 'FontName', 'Times New Roman');
xlabel('Frequency (Hz)','FontSize',20,'FontName','Times New Roman');
ylabel('Phase Difference ','FontSize',20,'FontName','Times New Roman');
grid on;

% Lag time diagram
subplot(3, 1, 2);
plot(frequencies, lag_time, 'LineWidth', 1.5);
set(gca,'FontSize',20, 'FontName', 'Times New Roman');
xlabel('Frequency (Hz)','FontSize',20,'FontName','Times New Roman');
ylabel('Lag Time ','FontSize',20,'FontName','Times New Roman');
grid on;

% Coherence diagram
subplot(3, 1, 3);
plot(f, Cxy, 'LineWidth', 1.5);
set(gca,'FontSize',20, 'FontName', 'Times New Roman');
xlabel('Frequency (Hz)','FontSize',20,'FontName','Times New Roman');
ylabel('Coherence','FontSize',20,'FontName','Times New Roman');
grid on;
