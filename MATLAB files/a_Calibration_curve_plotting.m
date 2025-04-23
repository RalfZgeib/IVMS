% Calibration curve
% Input measured voltage with calibration weights
g = 9.81;
m_calibration = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 150, 200, 250, 300, 400, 500]; % Known weights in g
F_calibration = m_calibration/1000 * g;

V_0=[15.4, 25.6, 33.3, 43.5, 53.5, 62, 71.3, 80.1, 89.7, 98.6, 141.3, 190.9, 236.2, 283.2, 386.6, 480.6];
V_1=[15.7, 26.1, 34, 44.4, 54.6, 63.2, 72.7, 81.7, 91.5, 100.6, 144.1, 194.7, 240.9, 288.9, 394.3, 490.2];
V_2=[16.3, 25.2, 34.2, 43.2, 52.4, 61.5, 72.7, 80.1, 89.8, 98.9, 145, 191.6, 245.2, 285.5, 388.2, 493.2];
V_3=[16.8, 25.4, 36.5, 46.6, 54.4, 63.7, 72.8, 81.9, 92.5, 100.4, 149.3, 197.5, 248.5, 294.1, 388.6, 501.7];

% Initial Voltage values
V_1_i=6.7;
V_2_i=6.5;
V_3_i=6.4;

% Tared Voltage values
V_0_T=V_0-V_1_i;
V_1_T=V_1-V_1_i;
V_2_T=V_2-V_2_i;
V_3_T=V_3-V_3_i;

% Calculating avg and std for V and F_measured
V=[V_1_T; V_2_T; V_3_T];
V_avg=mean(V,1);
V_std=std(V,0, 1);
F_std=V_std/1000*g;

% Calculate Linear Fit (Calibration Equation)
p = polyfit(V_avg, F_calibration, 1); % Fit line: Force = p(1)*Voltage + p(2)
F_fitted = polyval(p, V_avg); % Fitted values

% Calculate Error Envelope
residuals = F_calibration - F_fitted; % Residuals (difference between actual and fitted forces)
std_residuals = std(residuals); % Standard deviation of residuals
upper_bound = F_fitted + std_residuals; % Upper error bound
lower_bound = F_fitted - std_residuals; % Lower error bound

% Plotting Calibration Curves
figure;
hold on;

% Plot Data Points using Scatter (Older MATLAB Syntax)
h_data = scatter(V_avg, F_calibration, 16, ... % Size of the markers
    'MarkerEdgeColor', 'b', ... % Blue edge for points
    'MarkerFaceColor', 'b', ... % Blue fill for points
    'DisplayName', 'Measured Data');

% Plot Linear Fit
h_fit = plot(V_avg, F_fitted, '-', 'LineWidth', 1, 'Color', 'r', ...
    'DisplayName', 'Linear Model Fit');

% Plot Error Envelope
h_fill = fill([V_avg, fliplr(V_avg)], [upper_bound, fliplr(lower_bound)], [0.1, 0.1, 0.1], ...
    'EdgeColor', 'none', 'FaceAlpha', 0.2, 'DisplayName', 'STD Residual Error Envelope');

% Annotate Plot
xlabel('Measured Voltage (mV)');
ylabel('Calibration Force (N)');
title('Load Cell Calibration Curve');

% Custom Legend Order
legend([h_data, h_fit, h_fill], ...
    {'Measured Data', 'Linear Model Fit', 'STD Residual Error Envelope'}, ...
    'Location', 'Best');

% Calculate Average Percent Error
percent_errors = abs((F_calibration - F_fitted) ./ F_calibration) * 100;
avg_percent_error = mean(percent_errors);

% Display Calibration Equation, R^2, and Average Percent Error on the Graph
% Display the linear fit equation (p(1) is the slope, p(2) is the intercept)
annotation_text = sprintf('Linear Model Fit: F = %.5f * V + %.5f\nR^2 = %.5f\nAvg. Percent Error = %.2f%%', p(1), p(2), R_squared, avg_percent_error);
text(min(V_avg) + 0.1 * range(V_avg), max(F_calibration) - 0.1 * range(F_calibration), annotation_text, ...
    'FontSize', 10, 'BackgroundColor', 'white', 'EdgeColor', 'black', 'Margin', 5);

hold off;