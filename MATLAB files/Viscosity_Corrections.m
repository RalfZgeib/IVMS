classdef Viscosity_Corrections
    properties
        R % Capillary radii (m) - 1x3 array
        L % Capillary lengths (m) - 1x3 array
        Q % Flow rates (m³/s) - 6x1 array
        P % 6x3 matrix of total pressures (rows = flow rates, cols = needles)
        P_std % Total STD of measured pressures
        P_water % Total water pressures for friction correction
        P_water_STD % Total water pressure STDs for STD friction correction
        CY_params = [] % Carreau-Yasuda model parameters [eta_0, eta_inf, lambda, n, a]
        PL_params = [] % Power Law model parameters [K, n]
        initial_PL % Initial guess parameters for Power Law fitting
        initial_CY % Initial guess parameters for Carreau-Yasuda fitting
        do_FC % Flag whether to do Friction Correction or not
        do_BC % Flag whether to do Bagley Correction or not
        do_SC % Flag whether to do Slip Correction or not
        do_RC % Flag whether to do Robinowitsch Correction or not
    end

    methods
        % Constructor
        function obj = Viscosity_Corrections(R, L, Q, P, P_std, P_water, P_water_STD, initial_PL, initial_CY, do_FC, do_BC, do_SC, do_RC)
            obj.R = R;
            obj.L = L;
            obj.Q = Q;
            obj.P = P;
            obj.P_std = P_std;
            obj.P_water = P_water;
            obj.P_water_STD = P_water_STD;
            obj.initial_PL = initial_PL;
            obj.initial_CY = initial_CY;
            obj.do_FC = do_FC;
            obj.do_BC = do_BC;
            obj.do_SC = do_SC;
            obj.do_RC = do_RC;
        end

        function [P_F_corrected, P_friction_std] = apply_friction_correction(obj)
            % Input:
            % P_bioink: Matrix of bioink pressures (rows = flow rates, cols = needles)
            % P_water: Matrix of water pressures (rows = flow rates, cols = needles)
            if obj.do_FC == 1
                % Water viscosity at 25°C (Pa·s)
                eta_water = 0.891e-3;
                if any(size(obj.P) > 1)
                    % Calculate theoretical water pressure for each needle and flow rate
                    P_theoretical = zeros(size(obj.P_water));
                    for i = 1:length(obj.Q)
                        for j = 1:size(obj.P, 2)
                            gamma_dot = (4 * obj.Q(i)) / (pi * obj.R(j)^3);
                            tau_w = eta_water * gamma_dot;
                            P_theoretical(i,j) = (2 * obj.L(j) * tau_w) / obj.R(j);
                        end
                    end
                    % Calculate friction pressure
                    P_friction = obj.P_water - P_theoretical;
                    P_friction_std = obj.P_water_STD;

                else % Use needle 2 values for calculations if only 1 pressure column is given
                    gamma_dot = (4 * obj.Q) / (pi * obj.R(2)^3);
                    tau_w = eta_water * gamma_dot;
                    P_theoretical = (2 * obj.L(2) * tau_w) / obj.R(2);

                    % Calculate friction pressure
                    P_friction = obj.P_water(:,1) - P_theoretical;
                    P_friction_std = obj.P_water_STD(:,1);
                end
                
                % Subtract friction pressure from bioink pressures
                P_F_corrected = obj.P - P_friction;
                
                % Ensure no negative pressures
                P_F_corrected = max(P_F_corrected, 0);
           else
                P_F_corrected = obj.P; % Assuming no friction correction, pressures remain the same
                P_friction_std = obj.P_water_STD(:,1);
           end
        end


        % Bagley correction to calculate flow pressure (using R(2) and R(3))
        function [L_D, P_flow] = apply_bagley_correction(obj, P_F_corrected)
            % L/D for needles 2 and 3 of the same D
            L_D = obj.L(2:3) ./ (2 * obj.R(2));
            
            if any(size(obj.P) == 1)
                % Skip Bagley and slip corrections for single needle
                P_flow = P_F_corrected;  % Use direct pressure measurements
            else
                if obj.do_BC == 1
                    P_flow = zeros(size(obj.Q)); % Initialize flow pressure array for 19G/1.5" needle
                    % Fit pressure data for the current flow rate across needles 2 and 3
                    for i = 1:length(obj.Q)
                        Bagley_coeffs = polyfit(L_D, P_F_corrected(i, 2:3), 1); % Fit linear model
                        
                        % Compute P_flow = P_total - P_entrance
                        P_flow(i) = P_F_corrected(i, 2) - Bagley_coeffs(2); % Flow pressure array for 19G/1.5" needle only
                        
                        % Check if P_flow(i) is negative and report an error
                        if P_flow(i) < 0
                            fprintf(['P_flow at index %d is negative. This may indicate an issue with the input data...' ...
                                ' such as measured pressure with the longer L needle are smaller than that of the shorter L needle .\n'], i);
                            error('Negative P_flow encountered at index %d. Value: %f', i, P_flow(i));
                        end
                    end
                else
                    P_flow = P_F_corrected(:,2);
                end
            end
        end


        % Calculate apparent shear rate and shear stress
        function [gamma_dot_app, tau_w_app, eta_app] = calculate_apparent_params(obj, P_flow, Q)
            if any(size(obj.P) == 1)
                % Single needle calculations
                gamma_dot_app = (4 * Q) / (pi * obj.R(2)^3);
                tau_w_app = P_flow * obj.R(2) / (2 * obj.L(2));
                eta_app = tau_w_app ./ gamma_dot_app;
            else
            
            % Initialize & calculate gamma for each capillary
            gamma_dot_app = zeros(length(Q), length(obj.R));
            for i = 1:length(obj.R)
                gamma_dot_app(:,i) = (4 * Q) / (pi * obj.R(i)^3);  % Apparent shear rate (1/s)
            end 

            % Calculate parameters for capillary 2 (19G/1.5" needle)
            tau_w_app = P_flow * obj.R(2) / (2 * obj.L(2));  % Apparent wall shear stress (Pa)
            eta_app = tau_w_app ./ gamma_dot_app(:,2);  % Apparent viscosity (Pa·s)
            end
        end

        
        % Perform slip correction using 19G/22G needles of the same length but different radii
        function Q_corrected = apply_slip_correction(obj, gamma_dot_app)
            if obj.do_SC == 1
                if any(size(obj.P) > 1)
                    slip_data = zeros(length(obj.Q), 2); % Only use first two radii (same length)
                    slip_velocity = zeros(length(obj.Q),1); % Only calculated for capillary 2
                    Q_corrected = zeros(length(obj.Q),1); % Corrected flow for capillary 2
    
                    for r = 1:2 % radii of capillary 1 & 2
                        slip_data(:, r) = 1 ./ gamma_dot_app(:,r); % Inverse shear rate
                    end
    
                    % Fit Mooney's equation
                    for i = 1:length(obj.Q)
                        Mooney_coeffs = polyfit(1 ./ obj.R(1:2), slip_data(i,:), 1);
                        % Slip velocity = slope/2
                        slip_velocity(i) = Mooney_coeffs(1) / 2;
                        % Adjust flow rate for slip
                        Q_corrected(i) = obj.Q(i) - slip_velocity(i) * 2 * pi * obj.R(2)^2; % Corrected flow for capillary 2
                    end
                else
                    Q_corrected = obj.Q; % No correction if insufficient radii data
                end
            else
                Q_corrected = obj.Q;
            end
        end


        function [n, gamma_dot_true, eta_corrected] = rabinowitsch_correction(obj, gamma_dot_app, tau_w_app)

            if any(size(obj.P) == 1)
                % Use the single needle data directly
                n = 1;
                gamma_dot_true = gamma_dot_app;
                eta_corrected = tau_w_app ./ gamma_dot_app;
            else
                if obj.do_RC == 1
                    % Extract data for needle 2 only
                    gamma_dot = gamma_dot_app(:,2);
                    % Use log values to compute single n value across all data
                    log_gamma = log10(gamma_dot);
                    log_tau = log10(tau_w_app);
                    
                    % Fit a line to log-log data to get single power-law index
                    coeffs = polyfit(log_gamma, log_tau, 1);
                    n = coeffs(1);  % Slope of log-log plot is the power-law index
                    
                    % Compute true shear rate using single n value
                    gamma_dot_true = gamma_dot .* (3 * n + 1) ./ (4 * n);
                    eta_corrected = tau_w_app ./ gamma_dot_true;
                else
                    n = 1;
                    gamma_dot_true = gamma_dot_app(:,2);
                    eta_corrected = tau_w_app ./ gamma_dot_app(:,2);
                end
            end
        end


        function obj = fit_carreau_yasuda(obj, gamma_dot_true, eta_corrected)
            % Objective function for CY model fitting
            CY_model = @(params, gamma_dot) params(2) + (params(1) - params(2)) .* ...
                (1 + (params(3) .* gamma_dot).^params(4)).^((params(5) - 1) / params(4));
        
            % Parameter bounds [eta_0, eta_inf, lambda, a, n]
            lower_bounds = [0, 0, 0, 0, 0];
            upper_bounds = [inf, inf, inf, inf, 1];
        
            % Global optimization using Simulated Annealing
            fun = @(params) sum((eta_corrected - CY_model(params, gamma_dot_true)).^2);
        
            % Set options for simulated annealing
            sa_options = optimoptions('simulannealbnd', 'MaxIterations', 1e8, ...
                'InitialTemperature', 1.0, 'Display', 'off');
        
            % Perform simulated annealing
            [global_params, ~] = simulannealbnd(fun, obj.initial_CY, lower_bounds, upper_bounds, sa_options);
        
            % Use global optimization result as the initial guess for lsqcurvefit
            obj.initial_CY = global_params;
        
            % Set optimization options for lsqcurvefit
            options = optimoptions('lsqcurvefit', 'TolFun', 1e-8, 'MaxIterations', 1e8, ...
                'MaxFunctionEvaluations', 1e8, 'Display', 'off', 'Algorithm', 'trust-region-reflective');
        
            % Fit the model using lsqcurvefit
            fitted_params = lsqcurvefit(CY_model, obj.initial_CY, gamma_dot_true, eta_corrected, ...
                lower_bounds, upper_bounds, options);
        
            % Calculate adjusted R² value
            eta_pred = CY_model(fitted_params, gamma_dot_true);
            SS_res = sum((eta_corrected - eta_pred).^2);
            SS_tot = sum((eta_corrected - mean(eta_corrected)).^2);
            n_points = length(eta_corrected);
            n_params = 5;  % Number of CY model parameters
            df_res = n_points - n_params;  % Residual degrees of freedom
            df_tot = n_points - 1;         % Total degrees of freedom
            
            if df_res > 0 && df_tot > 0
                R2_adj = 1 - (SS_res / df_res) / (SS_tot / df_tot);
                R2_adj = max(0, min(1, R2_adj));  % Bound between 0 and 1
            else
                R2_adj = NaN;
            end
        
            % Store fitted parameters and R²
            obj.CY_params = [fitted_params, R2_adj]; % 5 parameters + 1 R² value
        end

        
        
        function obj = fit_power_law(obj, gamma_dot_true, eta_corrected)
            % Objective function for Power Law model fitting
            PL_model = @(params, gamma_dot) params(1) .* gamma_dot.^(params(2)-1);
        
            % Parameter bounds [K, n]
            lower_bounds = [eps, 0];
            upper_bounds = [1e6, 1];
        
            % Set optimization options
            options = optimoptions('lsqcurvefit', 'TolFun', 1e-8, 'MaxIterations', 10e8, ...
                'MaxFunctionEvaluations', 10e8, 'Display', 'off', 'Algorithm', 'trust-region-reflective', 'FiniteDifferenceStepSize', 1e-8);
        
            % Fit the model
            fitted_params = lsqcurvefit(PL_model, obj.initial_PL, gamma_dot_true, eta_corrected, ...
                lower_bounds, upper_bounds, options);
            
            % Calculate adjusted R² value
            eta_pred = PL_model(fitted_params, gamma_dot_true);
            SS_res = sum((eta_corrected - eta_pred).^2);
            SS_tot = sum((eta_corrected - mean(eta_corrected)).^2);
            n_points = length(eta_corrected);
            n_params = 2;  % Number of PL model parameters
            df_res = n_points - n_params;  % Residual degrees of freedom
            df_tot = n_points - 1;         % Total degrees of freedom
            
            if df_res > 0 && df_tot > 0
                R2_adj = 1 - (SS_res/df_res)/(SS_tot/df_tot);
                R2_adj = max(0, min(1, R2_adj));  % Bound between 0 and 1
            else
                R2_adj = NaN;
            end
    
            % Store fitted parameters and R²
            obj.PL_params = [fitted_params, R2_adj]; % 2 parameters + 1 R² value
        end


        % Calculate viscosity using Carreau-Yasuda model
        function eta_CY = carreau_yasuda_viscosity(obj, gamma_dot_true)
            % Ensure parameters are fitted
            if isempty(obj.CY_params)
                error('CY parameters not fitted. Run fit_carreau_yasuda first.');
            end
            params = obj.CY_params;
            eta_CY = params(2) + (params(1) - params(2)) .* ...
                (1 + (params(3) .* gamma_dot_true).^params(4)).^((params(5) - 1) / params(4));
        end


        function eta_PL = power_law_viscosity(obj, gamma_dot_true)
            % Ensure parameters are fitted
            if isempty(obj.PL_params)
                error('Power Law parameters not fitted. Run fit_power_law first.');
            end
            K = obj.PL_params(1);
            n = obj.PL_params(2);
            eta_PL = K * gamma_dot_true.^(n-1);
        end

        function [P_std_FC, eta_app_STD, eta_app_slip_STD, eta_corrected_STD] = compute_eta_corrected_std(obj, P_friction_std, gamma_dot_app, gamma_dot_app_slip, gamma_dot_true)
        % FRICTION_THEN_RABINOWITSCH_STD
            % Applies only friction correction, computes STD for the 
            % friction-corrected pressures, calculates the apparent viscosity STD, 
            % and then does a single Rabinowitsch correction. Returns only the 
            % final eta_corrected_STD.
        
            % 1) Apply friction correction
            if abs(obj.P_std - P_friction_std) > obj.P_std
                P_std_FC = obj.P_std;
            else
                P_std_FC = abs(obj.P_std - P_friction_std);
            end

            % 2) Compute "apparent" shear stress & viscosity and their STDs 
            %     (assuming Q is known exactly => gamma_dot_std = 0)
        
            % Apparent wall shear stress STD
            tau_w_app_STD = (P_std_FC(:,2) .* obj.R(2)) ./ (2 .* obj.L(2)); % (Pa)
        
            % Apparent viscosity STD
            eta_app_STD = tau_w_app_STD ./ gamma_dot_app(:,2); % (Pa·s)

            % Apparent slip viscosity STD
            eta_app_slip_STD = tau_w_app_STD ./ gamma_dot_app_slip(:,2); % (Pa·s)
        
            % 3) Apply Rabinowitsch correction same n-value calculated for eta_corrected
            eta_corrected_STD = tau_w_app_STD ./ gamma_dot_true; % (Pa·s)
        end


        function [results, obj] = analyze(obj)
            
            % Apply Friction Correction
            [P_F_corrected, P_friction_std] = obj.apply_friction_correction();

            % Apply Bagley correction
            [L_D, P_flow] = obj.apply_bagley_correction(P_F_corrected);
            
            % Calculate apparent parameters
            [gamma_dot_app, tau_w_app, eta_app] = obj.calculate_apparent_params(P_flow, obj.Q);
            
            % Apply slip correction
            Q_corrected = obj.apply_slip_correction(gamma_dot_app);
            
            % Recalculate parameters with slip-corrected flow rate
            [gamma_dot_app_slip, tau_w_app_slip, eta_app_slip] = obj.calculate_apparent_params(P_flow, Q_corrected);
            
            % Apply Rabinowitsch correction
            [n, gamma_dot_true, eta_corrected] = obj.rabinowitsch_correction(gamma_dot_app_slip, tau_w_app_slip);
            
            % Fit rheological models
            obj = obj.fit_power_law(gamma_dot_true, eta_corrected);
            obj = obj.fit_carreau_yasuda(gamma_dot_true, eta_corrected);
            
            % Calculate STDs
            [P_std_FC, eta_app_STD, eta_app_slip_STD, eta_corrected_STD] = compute_eta_corrected_std(obj, P_friction_std, gamma_dot_app, gamma_dot_app_slip, gamma_dot_true);

            % Store all results in a structure
            results.P_F_corrected = P_F_corrected;
            results.P_friction_std = P_friction_std;
            results.L_D = L_D;
            results.P_flow = P_flow;
            results.gamma_dot_app = gamma_dot_app;
            results.tau_w_app = tau_w_app;
            results.eta_app = eta_app;
            results.Q_corrected = Q_corrected;
            results.gamma_dot_app_slip = gamma_dot_app_slip;
            results.tau_w_app_slip = tau_w_app_slip;
            results.eta_app_slip = eta_app_slip;
            results.n = n;
            results.gamma_dot_true = gamma_dot_true;
            results.eta_corrected = eta_corrected;
            results.PL_params = obj.PL_params;
            results.CY_params = obj.CY_params;
            results.eta_app_STD = eta_app_STD;
            results.eta_app_slip_STD = eta_app_slip_STD;
            results.eta_corrected_STD = eta_corrected_STD;
            results.P_std_FC = P_std_FC;
        end


        function plot_corrections(obj, results)
            % Create figure with consistent formatting
            figure('Position', [100 100 1200 800]);
            t = tiledlayout(3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
            
            % Color scheme
            %colors = {'#0072BD', '#D95319', '#77AC30', '#800080'}; % MATLAB blue, orange, green, and purple
            colors = {[0, 0.4470, 0.7410], [0.8500, 0.3250, 0.0980], [0.4660, 0.6740, 0.1880], [0.4940, 0.1840, 0.5560]};
            shadesBlue = parula(length(obj.Q));  % Color-blind colormap
            markers = {'o', 's', 'd', '^', 'v', '>'};  % Six different markers for six flow rates
            
            % Plot 1: Original vs Friction Corrected Pressure for needle 2 
            nexttile
            hold on
            % Raw pressure data
            plot(obj.Q, obj.P(:,2), [markers{1} '-'], ...
                'Color', colors{4}, 'LineWidth', 1.5, ...
                'DisplayName', sprintf('Original P '));
            % STD
            errorenvelope(obj.Q, obj.P(:,2), obj.P_std(:,2), colors{4}, 0.1)
            
            % Friction corrected pressures
            plot(obj.Q, results.P_F_corrected(:,2), [markers{2} '--'], ...
                'Color', colors{4}, 'LineWidth', 1.5, ...
                'DisplayName', 'Friction Corrected P');
            % STD
            errorenvelope(obj.Q, results.P_F_corrected(:,2), results.P_std_FC(:,2), colors{4}, 0.1)
            
            xlabel('Flow Rate (m³/s)')
            ylabel('Pressure (Pa)')
            title('Original vs Friction Corrected Pressure')
            legend('Location', 'northwest', 'Box', 'off')
            grid on
            box on
            hold off

            % Plot 2: P vs. L/D for each flow rate Q
            nexttile
            hold on
            for i = 1:length(obj.Q)
                plot(results.L_D, results.P_F_corrected(i,2:3), [markers{i} '-'], ...
                    'Color', shadesBlue(i,:), 'LineWidth', 1.5, ...
                    'DisplayName', sprintf('Q = %.2e m³/s', obj.Q(i)));
            end
            xlabel('L/D Ratio')
            ylabel('Pressure Drop (Pa)')
            title({'Friction Corrected P vs. L/D', '(using needles with same D and different L for Bagley Correction)'})
            legend('Location', 'northwest', 'Box', 'off')
            grid on
            box on
            hold off
            
            % Plot 3: Original vs Bagley Corrected Pressure for needle 2
            nexttile
            hold on

            % Friction corrected pressures
            plot(obj.Q, results.P_F_corrected(:,2), [markers{1} '-'], ...
                'Color', colors{1}, 'LineWidth', 1.5, ...
                'DisplayName', sprintf('Original P (L/D = %.1f)', results.L_D(2)));
            % STD
            errorenvelope(obj.Q, results.P_F_corrected(:,2), results.P_std_FC(:,2), colors{1}, 0.1)

            % Bagley corrected pressures
            plot(obj.Q, results.P_flow, [markers{2} '--'], ...
                'Color', colors{1}, 'LineWidth', 1.5, ...
                'DisplayName', 'Bagley Corrected P');
            % STD
            errorenvelope(obj.Q, results.P_flow, results.P_std_FC(:,2), colors{1}, 0.1)

            xlabel('Flow Rate (m³/s)')
            ylabel('Pressure (Pa)')
            title('Original vs Bagely Corrected Pressure')
            legend('Location', 'northwest', 'Box', 'off')
            grid on
            box on
            hold off

            % Plot 4: Wall Slip Correction for needle 2
            nexttile
            hold on
            plot(obj.Q, results.gamma_dot_app(:,2), [markers{1} '-'], ...
                'Color', colors{2}, 'LineWidth', 1.5, ...
                'DisplayName', sprintf('Apparent γ̇ '));
            
            plot(results.Q_corrected, results.gamma_dot_app_slip(:,2), [markers{2} '--'], ...
                'Color', colors{2}, 'LineWidth', 1.5, ...
                'DisplayName', sprintf('Slip-corrected γ̇ '));
        
            xlabel('Flow Rate (m³/s)')
            ylabel('Shear Rate (s^{-1})')
            title('Wall Slip Correction using 19G/22G needles with the same L')
            legend('Location', 'northwest', 'Box', 'off')
            grid on
            box on
            hold off
            
            % Plot 5: Rabinowitsch Correction
            nexttile
            hold on
            
            % Slip-corrected viscosity
            plot(results.gamma_dot_app_slip(:,2), results.eta_app_slip, [markers{1} '-'], ...
                'Color', colors{3}, 'LineWidth', 1.5, ...
                'DisplayName', sprintf('Apparent Shear Rate & Viscosity'));
            % STD
            errorenvelope(results.gamma_dot_app_slip(:,2), results.eta_app_slip, results.eta_app_slip_STD, colors{3}, 0.1)

            % Robinowitsch corrected viscosity
            plot(results.gamma_dot_true, results.eta_corrected, [markers{2} '--'], ...
                'Color', colors{3}, 'LineWidth', 1.5, ...
                'DisplayName', sprintf('True Shear Rate & Viscosity'));
            % STD
            errorenvelope(results.gamma_dot_true, results.eta_corrected, results.eta_corrected_STD, colors{3}, 0.1)

            set(gca, 'XScale', 'log', 'YScale', 'log')
            xlabel('Shear Rate (s^{-1})')
            ylabel('Viscosity (Pa·s)')
            title('Rabinowitsch Correction')
            legend('Location', 'northeast', 'Box', 'off')
            grid on
            box on
            hold off
        
            % Plot 6: Plotting Rheology Model fitting
            nexttile
            hold on

            % Setting new colors for clarity
            colors_rheo = {'#77AC30', '#D95319', [0.4940, 0.1840, 0.5560]};
            
            % Plot viscosity with model fitting
            loglog(results.gamma_dot_true, results.eta_corrected, ...
                ['s' '--'], 'Color', colors_rheo{1}, 'LineWidth', 2, ...
                'DisplayName', sprintf('IVMS Measured Viscosity'));
            % STD
            errorenvelope(results.gamma_dot_true, results.eta_corrected, results.eta_corrected_STD, colors{3}, 0.1)

            % Generate smooth curves for models
            gamma_range = logspace(min(log10(results.gamma_dot_true)), ...
                max(log10(results.gamma_dot_true)), 100);
            
            % Plot Carreau-Yasuda fit
            eta_CY = obj.carreau_yasuda_viscosity(gamma_range);
            loglog(gamma_range, eta_CY, '-.', 'Color', colors_rheo{2}, 'LineWidth', 1.5, ...
                'DisplayName', sprintf('CY Model Fit (R² = %.5f)', obj.CY_params(6)));
            
            % Plot Power Law fit
            eta_PL = obj.power_law_viscosity(gamma_range);
            loglog(gamma_range, eta_PL, ':', 'Color', colors_rheo{3}, 'LineWidth', 1.5, ...
                'DisplayName', sprintf('PL Model Fit (R² = %.5f)', obj.PL_params(3)));
            
            xlabel('True Shear Rate (s^{-1})')
            ylabel('Corrected Viscosity (Pa·s)')
            title('Corrected IVMS Viscosity of 5% GelMA @ 25°C');
            legend('Location', 'northeast', 'Box', 'off')
            grid on
            box on
            set(gca, 'XScale', 'log', 'YScale', 'log');
            hold off
            
            % Overall formatting
            set(gcf, 'Color', 'white')
            set(findall(gcf,'-property','FontSize'),'FontSize', 12)
            set(findall(gcf,'-property','FontName'),'FontName', 'Arial')
        end
        
        function plot_rheology(obj, results)
            figure('Position', [100 100 800 600]);
            colors = {'#77AC30', '#D95319', [0.4940, 0.1840, 0.5560]};
            
            % Plot experimental data
            loglog(results.gamma_dot_true, results.eta_corrected, ...
                's', 'Color', colors{1}, 'LineWidth', 2, ...
                'DisplayName', sprintf('IVMS Measured Viscosity ±1σ'));
            hold on
            
            % Generate smooth curves for models
            gamma_range = logspace(min(log10(results.gamma_dot_true)), ...
                max(log10(results.gamma_dot_true)), 100);
            
            % Plot Carreau-Yasuda fit
            eta_CY = obj.carreau_yasuda_viscosity(gamma_range);
            loglog(gamma_range, eta_CY, '-.', 'Color', colors{2}, 'LineWidth', 1.5, ...
                'DisplayName', sprintf('CY Model Fit (R² = %.5f)', obj.CY_params(6)));
            
            % Plot Power Law fit
            eta_PL = obj.power_law_viscosity(gamma_range);
            loglog(gamma_range, eta_PL, ':', 'Color', colors{3}, 'LineWidth', 1.5, ...
                'DisplayName', sprintf('PL Model Fit (R² = %.5f)', obj.PL_params(3)));
            
            xlabel('True Shear Rate (s^{-1})')
            ylabel('Corrected Viscosity (Pa·s)')
            title({'Corrected IVMS Viscosity', '5% GelMA @ 25°C'});
            legend('Location', 'best', 'Box', 'off', 'FontSize', 20)
            grid on
            
            set(gcf, 'Color', 'white')
            set(gca, 'FontSize', 25, 'FontName', 'Arial')
        end
    
        
        function plot_comparisons(obj, results_array, labels, conditions)
            % Input parameters:
            % results_array: Cell array of results structures from different samples
            % labels: Cell array of sample labels (e.g., {'5%', '7%', '10%'})
            % conditions: Structure with fields:
            %   - concentrations: Cell array of concentrations to compare
            %   - temperatures: Cell array of temperatures to compare
            %
            % Creates a 1x2 layout comparing concentration effects (left) 
            % and temperature effects (right), each with an optional ±1σ envelope 
            % if 'eta_corrected_STD' is present in the results.

            % Create figure with 1x2 layout
            figure('Position', [100 100 1200 500]);
            t = tiledlayout(1, 2, 'TileSpacing', 'loose', 'Padding', 'compact');
            
            % Color schemes
            colors_conc = { [0, 0.4470, 0.7410], [0.9290, 0.6940, 0.1250], [0.4660, 0.6740, 0.1880] };
            colors_temp = { [0.9290, 0.6940, 0.1250], [0.8500, 0.3250, 0.0980], [0.4940, 0.1840, 0.5560] };
            markers = {'o', 's', 'd', '^', 'v'};
            
            % Plot 1: Concentration effect
            nexttile;
            hold on;
            for i = 1:length(conditions.concentrations)
                idx = strcmp(labels, conditions.concentrations{i});
                if any(idx)
                    % Extract data
                    gamma_vals = results_array{idx}.gamma_dot_true;
                    eta_vals = results_array{idx}.eta_corrected;
                    eta_std = results_array{idx}.eta_corrected_STD;
            
                    % Plot main line
                    loglog(gamma_vals, eta_vals, [markers{i} '-'], 'Color', colors_conc{i}, 'LineWidth', 1.5, ...
                        'DisplayName', [conditions.concentrations{i}, ' ±1σ']);
            
                    % Plot error envelope
                    errorenvelope(gamma_vals, eta_vals, eta_std, colors_conc{i}, 0.1);
                end
            end
            xlabel('True Shear Rate (s^{-1})');
            ylabel('Corrected Viscosity (Pa·s)');
            title('Concentration Effect');
            legend('Location', 'northeast', 'Box', 'off', 'FontSize', 20);
            grid on;
            box on;
            hold off; % Reset for next plot
            
            % Plot 2: Temperature effect
            nexttile;
            hold on;
            for i = 1:length(conditions.temperatures)
                idx = strcmp(labels, conditions.temperatures{i});
                if any(idx)
                    % Extract data
                    gamma_vals = results_array{idx}.gamma_dot_true;
                    eta_vals = results_array{idx}.eta_corrected;
                    eta_std = results_array{idx}.eta_corrected_STD;
            
                    % Plot main line
                    loglog(gamma_vals, eta_vals, [markers{i} '-'], 'Color', colors_temp{i}, 'LineWidth', 1.5, ...
                        'DisplayName', [conditions.temperatures{i}, ' ±1σ']);
            
                    % Plot error envelope
                    errorenvelope(gamma_vals, eta_vals, eta_std, colors_temp{i}, 0.1);
                end
            end
            xlabel('True Shear Rate (s^{-1})');
            ylabel('Corrected Viscosity (Pa·s)');
            title('Temperature Effect');
            legend('Location', 'northeast', 'Box', 'off', 'FontSize', 20);
            grid on;
            box on;
            hold off;
            
            % Apply global formatting
            all_axes = findall(gcf, 'Type', 'axes');
            set(all_axes, 'XScale', 'log', 'YScale', 'log');
            set(gcf, 'Color', 'white');
            set(findall(gcf, '-property', 'FontSize'), 'FontSize', 25);
            set(findall(gcf, '-property', 'FontName'), 'FontName', 'Arial');
        end


        function plot_all_viscosities(obj, results_array, labels)
            % Create figure with professional formatting
            figure('Position', [100 100 800 600]);
            
            % Define color scheme and markers
            % Define custom color scheme
            colors = {
                [0, 0.4470, 0.7410],     % MATLAB blue for 5% @25C
                [0.9290, 0.6940, 0.1250], % MATLAB yellow for 7% @25C
                [0.8500, 0.3250, 0.0980], % MATLAB orange for 7% @30C
                [0.4940, 0.1840, 0.5560], % MATLAB purple for 7% @37C
                [0.4660, 0.6740, 0.1880]  % MATLAB green for 10% @25C
            };

            %colors = turbo(length(results_array));
            markers = {'o', 's', 'd', '^', 'v'};
            
            hold on
            for i = 1:length(results_array)
                gamma_vals = results_array{i}.gamma_dot_true;
                eta_vals = results_array{i}.eta_corrected;
                eta_std = results_array{i}.eta_corrected_STD;
       
                % Plot the mean curve
                loglog(gamma_vals, eta_vals, ...
                    [markers{i} '-'], ...
                    'Color', colors{i}, ...
                    'LineWidth', 1.5, ...
                    'DisplayName', [labels{i}, ' ±1σ']);
                
                % If the standard deviations exist, plot ±1σ envelope
                if isfield(results_array{i}, 'eta_corrected_STD') && ~isempty(results_array{i}.eta_corrected_STD)
                    errorenvelope(gamma_vals, eta_vals, eta_std, colors{i}, 0.1)
                end
            end
            
            % Formatting
            xlabel('True Shear Rate (s^{-1})')
            ylabel('Corrected Viscosity (Pa·s)')
            title('Measured GelMA Viscosity')
            legend('Location', 'northeast', 'Box', 'off', 'FontSize', 20)
            grid on
            box on
            
            % Set axes properties
            set(gca, 'XScale', 'log', 'YScale', 'log')
            
            % Overall formatting
            set(gcf, 'Color', 'white')
            set(findall(gcf,'-property','FontSize'),'FontSize', 25)
            set(findall(gcf,'-property','FontName'),'FontName', 'Arial')
        end

        function plot_rheometer_comparison_std(obj, results, rheometer_file, sample_name)
            % Load and process rheometer data from CSV
            rheometer_raw = readmatrix(rheometer_file);
            rheometer_data.gamma_dot = rheometer_raw(:,1);
            rheometer_data.eta = rheometer_raw(:,2);
        
            % Calculate percent error at each capillary data point
            percent_errors = zeros(length(results.gamma_dot_true), 1);
            for i = 1:length(results.gamma_dot_true)
                percent_errors(i) = abs( (results.eta_corrected(i) - rheometer_data.eta(i)) ...
                                         / rheometer_data.eta(i) ) * 100;
            end
        
            mean_error = mean(percent_errors);
            std_error = std(percent_errors);
        
            % Create figure with professional formatting
            figure('Position', [100 100 800 600]);
            colors = {[119, 172, 48] / 255, [217, 83, 25] / 255, [0.4940, 0.1840, 0.5560], [0, 114, 189] / 255}; % green, orange, purple, blue

        
            % Plot experimental data
            gamma_vals = results.gamma_dot_true;
            eta_vals   = results.eta_corrected;
            eta_std = results.eta_corrected_STD;

            loglog(gamma_vals, eta_vals, 's', ...
                'Color', colors{1}, 'LineWidth', 2, ...
                'DisplayName', sprintf('IVMS Measurements ±1σ (Mean Error: %.2f%% ± %.2f%%)', ...
                mean_error, std_error));
            hold on
        
            % Add error envelope around experimental data
            errorenvelope(gamma_vals, eta_vals, eta_std, colors{1}, 0.1)
        
            % Plot rheometer data
            loglog(rheometer_data.gamma_dot, rheometer_data.eta, 'o', ...
                'Color', colors{4}, 'LineWidth', 2, ...
                'DisplayName', 'Rheometer Data');
        
            % Generate smooth curves for model fitting
            gamma_range = logspace(min(log10([results.gamma_dot_true; rheometer_data.gamma_dot])), ...
                max(log10([results.gamma_dot_true; rheometer_data.gamma_dot])), 100)';
        
            % --------------------------
            % Fit & plot: Carreau-Yasuda (Rheometer)
            % --------------------------
            obj = fit_carreau_yasuda(obj, rheometer_data.gamma_dot, rheometer_data.eta);
            eta_CY_rheo = obj.carreau_yasuda_viscosity(gamma_range);
            loglog(gamma_range, eta_CY_rheo, '--', 'Color', colors{2}, ...
                'LineWidth', 1.5, 'DisplayName', ...
                sprintf('CY Model - Rheometer (R^2 = %.5f)', obj.CY_params(6)));
        
            % --------------------------
            % Fit & plot: Power Law (Rheometer)
            % --------------------------
            obj = fit_power_law(obj, rheometer_data.gamma_dot, rheometer_data.eta);
            eta_PL_rheo = obj.power_law_viscosity(gamma_range);
            loglog(gamma_range, eta_PL_rheo, '--', 'Color', colors{3}, ...
                'LineWidth', 1.5, 'DisplayName', ...
                sprintf('PL Model - Rheometer (R^2 = %.5f)', obj.PL_params(3)));
        
            % --------------------------
            % Fit & plot: Carreau-Yasuda (Capillary)
            % --------------------------
            obj = fit_carreau_yasuda(obj, results.gamma_dot_true, results.eta_corrected);
            eta_CY_cap = obj.carreau_yasuda_viscosity(gamma_range);
            loglog(gamma_range, eta_CY_cap, '-.', 'Color', colors{2}, ...
                'LineWidth', 1.5, 'DisplayName', ...
                sprintf('CY Model - Capillary (R^2 = %.5f)', obj.CY_params(6)));
        
            % --------------------------
            % Fit & plot: Power Law (Capillary)
            % --------------------------
            obj = fit_power_law(obj, results.gamma_dot_true, results.eta_corrected);
            eta_PL_cap = obj.power_law_viscosity(gamma_range);
            loglog(gamma_range, eta_PL_cap, ':', 'Color', colors{3}, ...
                'LineWidth', 1.5, 'DisplayName', ...
                sprintf('PL Model - Capillary (R^2 = %.5f)', obj.PL_params(3)));
        
            % Formatting
            xlabel('Shear Rate (s^{-1})')
            ylabel('Viscosity (Pa·s)')
            title(sample_name)
            legend('Location', 'northeast', 'Box', 'off', 'FontSize', 20)
            grid on
            box on
            hold off
            
            % Overall formatting
            set(gca, 'XScale', 'log', 'YScale', 'log')
            set(gcf, 'Color', 'white')
            set(findall(gcf,'-property','FontSize'),'FontSize', 25)
            set(findall(gcf,'-property','FontName'),'FontName', 'Arial')
        end
    end
end