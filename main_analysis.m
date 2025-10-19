% =========================================================================
% TCABR MOMENTUM ANALYSIS - MAIN SCRIPT
% =========================================================================

clc; clear; close all;
fprintf('Initiating TCABR momentum analysis script...\n\n');

% Add function directories to the MATLAB path
addpath('src', 'plotting', 'utils');

% --- STAGE 1: Setup and Data Loading ---
constants = setup_constants();
constants.analysis.num_iterations = 5000; % Define number of iterations here
exp_data = load_experimental_data(constants);

% --- STAGE 2: Velocity Profile Analysis ---
velocity_results = analyze_velocity_profile(exp_data, constants);

% --- Verification Step ---
fprintf('Verification: Plotting key velocity results...\n');

figure('Name', 'Velocity Profile Analysis Verification');
hold on;
grid on;
box on;

% Plot experimental data points with error bars (converted back to km/s for plotting)
errorbar(exp_data.r_Vphi_exp / constants.machine.a, ...
         exp_data.Vphi_exp, ...
         exp_data.Vphi_err_exp, ...
         'ko', 'MarkerFaceColor', 'k', 'DisplayName', 'Experimental Data');

% Plot the mean polynomial fit (converted to km/s)
plot(velocity_results.r_fine / constants.machine.a, ...
     velocity_results.poly_fit_avg / 1000, ...
     'b--', 'LineWidth', 2, 'DisplayName', 'Mean Polynomial Fit');

% Plot the 95% confidence interval band (converted to km/s)
fill([velocity_results.r_fine / constants.machine.a; ...
      flipud(velocity_results.r_fine / constants.machine.a)], ...
     [velocity_results.poly_fit_ci_lower / 1000; ...
      flipud(velocity_results.poly_fit_ci_upper / 1000)], ...
     'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', '95% Confidence Interval');

% Plot the Fourier-Bessel series fit (converted to km/s)
plot(velocity_results.r_fine / constants.machine.a, ...
     velocity_results.bessel_fit / 1000, ...
     'r-', 'LineWidth', 2, 'DisplayName', 'Fourier-Bessel Fit');

xlabel('Normalised Radius (r/a)');
ylabel('Toroidal Velocity, V_{\phi} (km/s)');
title('Velocity Profile Analysis Results');
legend('show', 'Location', 'best');
hold off;

fprintf('Script finished.\n');