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

% --- STAGE 2: Profile Analyses ---
velocity_results = analyze_velocity_profile(exp_data, constants);
temperature_results = analyze_temperature_profile(exp_data, constants);

% --- STAGE 3: Physics Profile Calculations ---
% Using the fine radial grid from the velocity analysis results
r_fine = velocity_results.r_fine;
magnetic_field = compute_magnetic_field(r_fine, constants);


% --- Verification Step ---
fprintf('Verification: Plotting key analysis results...\n');

% Figure 1: Velocity Profile Verification
% ... (código do plot de velocidade) ...
figure('Name', 'Velocity Profile Analysis Verification');
hold on; grid on; box on;
errorbar(exp_data.r_Vphi_exp / constants.machine.a, exp_data.Vphi_exp, exp_data.Vphi_err_exp, 'ko', 'MarkerFaceColor', 'k', 'DisplayName', 'Experimental Data');
plot(velocity_results.r_fine / constants.machine.a, velocity_results.poly_fit_avg / 1000, 'b--', 'LineWidth', 2, 'DisplayName', 'Mean Polynomial Fit');
fill([velocity_results.r_fine / constants.machine.a; flipud(velocity_results.r_fine / constants.machine.a)], [velocity_results.poly_fit_ci_lower / 1000; flipud(velocity_results.poly_fit_ci_upper / 1000)], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', '95% Confidence Interval');
plot(velocity_results.r_fine / constants.machine.a, velocity_results.bessel_fit / 1000, 'r-', 'LineWidth', 2, 'DisplayName', 'Fourier-Bessel Fit');
xlabel('Normalised Radius (r/a)'); ylabel('Toroidal Velocity, V_{\phi} (km/s)');
title('Velocity Profile Analysis Results'); legend('show', 'Location', 'best');
hold off;

% Figure 2: Ion Temperature Profile Verification
% ... (código do plot de temperatura) ...
figure('Name', 'Ion Temperature Analysis Verification');
hold on; grid on; box on;
errorbar(exp_data.r_Ti_exp / constants.machine.a, exp_data.Ti_reconstructed_exp, exp_data.Ti_reconstructed_err_exp, 'ko', 'MarkerFaceColor', 'k', 'DisplayName', 'Reconstructed Data');
plot(temperature_results.r_fine / constants.machine.a, temperature_results.profile_avg, 'r-', 'LineWidth', 2.5, 'DisplayName', 'Mean Canonical Fit');
fill([temperature_results.r_fine / constants.machine.a; flipud(temperature_results.r_fine / constants.machine.a)], [temperature_results.profile_ci_lower; flipud(temperature_results.profile_ci_upper)], 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', '95% Confidence Band');
xlabel('Normalised Radius (r/a)'); ylabel('Ion Temperature, T_{i} (eV)');
title('Ion Temperature Profile Analysis Results'); legend('show', 'Location', 'best'); xlim([0, 1]);
hold off;

% Figure 3: Magnetic Field Profile Verification
figure('Name', 'Magnetic Field Profile Verification');
hold on; grid on; box on;
plot(magnetic_field.r_fine / constants.machine.a, magnetic_field.total, 'k-', 'LineWidth', 2, 'DisplayName', 'B_{total}');
plot(magnetic_field.r_fine / constants.machine.a, magnetic_field.toroidal, 'b--', 'LineWidth', 2, 'DisplayName', 'B_{toroidal}');
plot(magnetic_field.r_fine / constants.machine.a, magnetic_field.poloidal, 'r-.', 'LineWidth', 2, 'DisplayName', 'B_{poloidal}');
xlabel('Normalised Radius (r/a)');
ylabel('Magnetic Field (T)');
title('Magnetic Field Profiles');
legend('show', 'Location', 'best');
hold off;

fprintf('Script finished.\n');