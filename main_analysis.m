% =========================================================================
% TCABR MOMENTUM ANALYSIS - MAIN SCRIPT
% =========================================================================

clc; clear; close all;
fprintf('Initiating TCABR momentum analysis script...\n\n');

% Add function directories to the MATLAB path
addpath('src', 'plotting', 'utils');

% --- STAGE 1: Setup and Data Loading ---
constants = setup_constants();
exp_data = load_experimental_data(constants);

% --- STAGE 2: Velocity Profile Analysis ---
velocity_results = analyze_velocity_profile(exp_data, constants);

% --- STAGE 3: Ion Temperature Profile Analysis ---
temperature_results = analyze_temperature_profile(exp_data, constants);


% --- Verification Step ---
fprintf('Verification: Plotting key analysis results...\n');

% Figure 1: Velocity Profile Verification
figure('Name', 'Velocity Profile Analysis Verification');
% ... (código do plot de velocidade da etapa anterior) ...
hold on; grid on; box on;
errorbar(exp_data.r_Vphi_exp / constants.machine.a, exp_data.Vphi_exp, exp_data.Vphi_err_exp, 'ko', 'MarkerFaceColor', 'k', 'DisplayName', 'Experimental Data');
plot(velocity_results.r_fine / constants.machine.a, velocity_results.poly_fit_avg / 1000, 'b--', 'LineWidth', 2, 'DisplayName', 'Mean Polynomial Fit');
fill([velocity_results.r_fine / constants.machine.a; flipud(velocity_results.r_fine / constants.machine.a)], [velocity_results.poly_fit_ci_lower / 1000; flipud(velocity_results.poly_fit_ci_upper / 1000)], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', '95% Confidence Interval');
plot(velocity_results.r_fine / constants.machine.a, velocity_results.bessel_fit / 1000, 'r-', 'LineWidth', 2, 'DisplayName', 'Fourier-Bessel Fit');
xlabel('Normalised Radius (r/a)'); ylabel('Toroidal Velocity, V_{\phi} (km/s)');
title('Velocity Profile Analysis Results'); legend('show', 'Location', 'best');
hold off;

% Figure 2: Ion Temperature Profile Verification
figure('Name', 'Ion Temperature Analysis Verification');
hold on; grid on; box on;

% Plot the "reconstructed" experimental data with its error bars
errorbar(exp_data.r_Ti_exp / constants.machine.a, ...
         exp_data.Ti_reconstructed_exp, ...
         exp_data.Ti_reconstructed_err_exp, ...
         'ko', 'MarkerFaceColor', 'k', 'DisplayName', 'Reconstructed Data');

% Plot the mean canonical fit
plot(temperature_results.r_fine / constants.machine.a, ...
     temperature_results.profile_avg, ...
     'r-', 'LineWidth', 2.5, 'DisplayName', 'Mean Canonical Fit');

% Plot the 95% confidence band
fill([temperature_results.r_fine / constants.machine.a; flipud(temperature_results.r_fine / constants.machine.a)], ...
     [temperature_results.profile_ci_lower; flipud(temperature_results.profile_ci_upper)], ...
     'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', '95% Confidence Band');

% Add text box with fit parameters
param_text = sprintf('Optimised Parameters:\nT_{i,0} = (%.1f \\pm %.1f) eV\nT_{i,a} = (%.1f \\pm %.1f) eV\n\\sigma = %.3f\nR^2 = %.3f', ...
    temperature_results.params_optimized.T0, temperature_results.param_uncertainties.T0_err, ...
    temperature_results.params_optimized.Ta, temperature_results.param_uncertainties.Ta_err, ...
    temperature_results.params_optimized.sigma, temperature_results.r_squared);
annotation('textbox', [0.2, 0.6, 0.3, 0.3], 'String', param_text, 'FitBoxToText', 'on', 'BackgroundColor', 'white', 'EdgeColor', 'k');

xlabel('Normalised Radius (r/a)');
ylabel('Ion Temperature, T_{i} (eV)');
title('Ion Temperature Profile Analysis Results');
legend('show', 'Location', 'best');
xlim([0, 1]);
hold off;

fprintf('Script finished.\n');