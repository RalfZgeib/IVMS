clear;
a = arduino('COM3');

% Initializations
t_final = 80; % Total duration (s)
Fs = 25; % Sampling frequency (Hz)
Ts = 1/Fs; % Sampling period (s)
Ns = t_final * Fs; % Total number of samples

% Calibration curve parameters F=p_F(1)⋅V+p_F(2) [V is in mV and F in N]
p_F=[0.01013, 0.02281];
g = 9.81; % Gravitational acceleration (m/s^2)

% Syringe Geometry
R_3mL=4.25; %3mL BD syringe barrel inner diameter
R_5mL=6; %5mL BD syringe barrel inner diameter

% Needle parameters (Comment/uncomment as needed)
R_Needle = 0.813/2; % Brown 19G needle inner radius in mm
%R_Needle = 0.508/2; % Black 22G needle inner radius in mm

% Preallocate arrays
P = zeros(Ns, 1); % Pressure array
P_cc = zeros(Ns, 1); % Pressure from calibration curve
time = zeros(Ns, 1); % Time array

% Set up plot for pressure
figure;
hPlotP = plot(time, P, 'LineWidth', 2.0);
xlabel('Time (s)');
ylabel('Pressure (Pa)');
grid on;

tic; % Start timing
for i = 1:Ns
    loopStartTime = toc;

    % Calculate pressure directly from voltage reading in Pa
    P(i) = (((readVoltage(a, 'A2') * 1000 - 2500) * 2) / 1000 * g) / (pi * R_3mL^2)*10^6;
    
    % Calculate pressure using calibration curve in Pa
    P_cc(i) = (p_F(1)*((readVoltage(a, 'A2') * 1000 - 2500) * 2)+p_F(2)) / (pi * R_3mL^2)*10^6;

    % Update time vector with actual elapsed time
    time(i) = toc;

    % Update plot data for pressure
    set(hPlotP, 'XData', time(1:i), 'YData', P(1:i));
    drawnow;

    % Calculate elapsed time for this iteration and adjust pause accordingly
    elapsedTime = toc - loopStartTime;
    pauseDuration = max(0, Ts - elapsedTime);
    pause(pauseDuration);
end

% Calculate timing & Acquisition rate
TotalelapsedTime = toc;
averageTimePerIteration = TotalelapsedTime / Ns;
acquisitionRate = 1 / averageTimePerIteration;

% Define headers for Excel output
headers = {'Time (s)', 'Pressure (Pa)', 'Calibration Curve Pressure (Pa)'};

% Define data matrix with only time and pressure
dataMatrix = [time P P_cc];

% Write headers to Excel
% SAMPLE name sample as measured
% writecell(headers, "SAMPLE_25C_19G_1.5in.xlsx", 'Range', 'A1');
% writecell(headers, "SAMPLE_25C_22G_1.5in.xlsx", 'Range', 'A1');
% writecell(headers, "SAMPLE_25C_22G_2in.xlsx", 'Range', 'A1');