% Clear workspace and command window
clear; 
clc; 
close all;
warning('off', 'MATLAB:ode45:IntegrationTolNotMet');
warning('off', 'MATLAB:ode23s:IntegrationTolNotMet');

% % Start parallel pool
% if isempty(gcp('nocreate'))
%     parpool('local');
% end

% Load glucose and heart rate data
glucose_data = readtable('PK0006_smoothed_glucose.csv');
heart_rate_data = readtable('PK0006_smoothed_hr.csv');

% Initial conditions
G0 = glucose_data.Average_Glucose(1);
H0 = heart_rate_data.hr_Average_Value(1);

% Baseline minimums for rest state
G_sleep = min(glucose_data.Average_Glucose);
H_sleep = min(heart_rate_data.hr_Average_Value);

% Time settings
t_eval = glucose_data.Time_Since_Start;
t_span = [0, max(t_eval)];

% Wake spike function
function spike = wake_spike(t, spike_amp, spike_width)
    tc = 260;
    if t >= tc
        spike = spike_amp * exp(-((t - tc)^2) / (2 * spike_width^2));
    else
        spike = 0;
    end
end

% ODE system
function dydt = glucose_heart_ODE(t, y, params, G_sleep, H_sleep)
    G = y(1);
    H = y(2);
    R_g = params(1);
    a_g = params(2);
    k_g = params(3);
    K_g = params(4);
    E_g = params(5);
    beta_H = params(6);
    H_v = params(7);
    spike_amp = params(8);
    spike_width = params(9);

    spike = wake_spike(t, spike_amp, spike_width);

    dG_dt = R_g + a_g*G*((H-H_sleep)/H_sleep) - k_g*(G/(G+K_g)) - E_g*G + spike;
    dH_dt = -beta_H*(H-H_sleep) + H_v*(G-G_sleep)/G_sleep;
    dydt = [dG_dt; dH_dt];
end

% Loss function
function total_loss = loss_function(params, glucose_data, heart_rate_data, G0, H0, t_eval, t_span, G_sleep, H_sleep)
    try
        opts = odeset('RelTol',1e-6,'AbsTol',1e-8,'MaxStep',1);
        sol = ode45(@(t,y) glucose_heart_ODE(t,y,params,G_sleep,H_sleep), t_span, [G0; H0], opts);

        % Check if the solver returned valid data
        if isempty(sol.x) || any(isnan(sol.y(:))) || any(isinf(sol.y(:))) || max(sol.x) < t_span(2)
            total_loss = 1e6; % Penalize incomplete solutions
            return;
        end

        G_model = interp1(sol.x, sol.y(1,:), t_eval, 'linear', 'extrap');
        H_model = interp1(sol.x, sol.y(2,:), t_eval, 'linear', 'extrap');

        % Check interpolated values too
        if any(isnan(G_model)) || any(isnan(H_model)) || any(isinf(G_model)) || any(isinf(H_model))
            total_loss = 1e6;
            return;
        end

        glucose_loss = mean((G_model - glucose_data.Average_Glucose).^2);
        heart_rate_loss = mean((H_model - heart_rate_data.hr_Average_Value).^2);
        total_loss = glucose_loss + heart_rate_loss;

        if isnan(total_loss) || isinf(total_loss)
            total_loss = 1e6;
        end
        
    catch
        % Hard solver crash
        total_loss = 1e6;
    end
end


% Optimization settings
param_bounds = [
    0.001, 10;    % R_g
    0.000, 0.1;    % a_g
    0.001, 10;    % k_g
    1,     10;    % K_g
    0.001, 5;     % E_g
    0.001, 10;     % beta_H
    0.001, 30;     % H_v
    0.01, 20;     % spike_amp
    10,   300     % spike_width
];
lb = param_bounds(:,1);
ub = param_bounds(:,2);

objective = @(params) loss_function(params, glucose_data, heart_rate_data, G0, H0, t_eval, t_span, G_sleep, H_sleep);

% GA (global optimization) settings
% options_ga = optimoptions('ga', ...
%     'Display', 'iter', ...
%     'MaxGenerations', 1000, ...
%     'PopulationSize', 60, ...
%     'UseParallel', false, ...
%     'FunctionTolerance', 1e-6);
% 
% [param_final, fval_final] = ga(objective, 9, [], [], [], [], lb, ub, [], options_ga);

% Local optimization settings
options_pso = optimoptions('particleswarm', ...
    'Display', 'iter', ...
    'UseParallel', false, ...
    'SwarmSize', 60, ...
    'MaxIterations', 300);

[param_final, fval_pso] = particleswarm(objective, 9, lb, ub, options_pso);

% % Local refinement
% local_options = optimoptions('fmincon', ...
%     'Algorithm', 'interior-point', ...
%     'SpecifyObjectiveGradient', false, ...
%     'Display', 'iter', ...
%     'MaxIterations', 2000, ...
%     'OptimalityTolerance', 1e-6, ...
%     'StepTolerance', 1e-8);
% 
% [param_final, fval_refined] = fmincon(objective, param_final, [], [], [], [], lb, ub, [], local_options);
% 

% Solve ODE with final parameters
opts = odeset('RelTol',1e-6,'AbsTol',1e-8);
sol_final = ode45(@(t,y) glucose_heart_ODE(t,y,param_final,G_sleep,H_sleep), t_span, [G0; H0], opts);

G_model = interp1(sol_final.x, sol_final.y(1,:), t_eval, 'linear', 'extrap');
H_model = interp1(sol_final.x, sol_final.y(2,:), t_eval, 'linear', 'extrap');

%% Plot results
figure;
plot(glucose_data.Time_Since_Start, glucose_data.Average_Glucose, 'b--', 'LineWidth', 2); hold on;
plot(t_eval, G_model, 'r-', 'LineWidth', 2);
xlabel('Time (min)', 'FontSize', 12);
ylabel('Glucose (mg/dl)', 'FontSize', 12);
title('Glucose', 'FontSize', 14);
legend('Observed', 'Fitted', 'Location', 'Best');

figure;
plot(heart_rate_data.Time_Since_Start, heart_rate_data.hr_Average_Value, 'g--', 'LineWidth', 2); hold on;
plot(t_eval, H_model, 'm-', 'LineWidth', 2);
xlabel('Time (min)', 'FontSize', 12);
ylabel('Heart Rate (BPM)', 'FontSize', 12);
title('Heart rate', 'FontSize', 14);
legend('Observed', 'Fitted', 'Location', 'Best');



% Print optimized parameters
fprintf('\nOptimized Parameters after Global + Local Optimization:\n');
fprintf('R_g         = %.3f\n', param_final(1));
fprintf('a_g         = %.3f\n', param_final(2));
fprintf('k_g         = %.3f\n', param_final(3));
fprintf('K_g         = %.3f\n', param_final(4));
fprintf('E_g         = %.3f\n', param_final(5));
fprintf('beta_H      = %.3f\n', param_final(6));
fprintf('H_v         = %.3f\n', param_final(7));
fprintf('spike_amp   = %.3f\n', param_final(8));
fprintf('spike_width = %.3f\n', param_final(9));


%%
% Sensitivity analysis
clc
DRIVER_sensitivity(param_final,G_sleep,H_sleep,G0,H0,t_span)