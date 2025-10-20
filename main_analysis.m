% =========================================================================
% TCABR MOMENTUM ANALYSIS - MAIN SCRIPT
% =========================================================================

clc; clear; close all;
fprintf('Initiating TCABR momentum analysis script...\n\n');

addpath('src', 'plotting', 'utils');

% --- STAGE 1: Setup and Data Loading ---
constants = setup_constants();
exp_data = load_experimental_data(constants);

r_fine = exp_data.r_profiles;

% --- STAGE 2: Profile Analyses ---
velocity_results = analyze_velocity_profile(exp_data, constants, r_fine);
temperature_results = analyze_temperature_profile(exp_data, constants, r_fine);

% --- STAGE 3: Physics Profile Calculations ---
magnetic_field = compute_magnetic_field(r_fine, constants);
derived_profiles = compute_derived_profiles(exp_data, temperature_results, magnetic_field, constants, r_fine);

% --- STAGE 4: Theoretical Model Calculations ---
theoretical_models = compute_theoretical_models(temperature_results, magnetic_field, derived_profiles, constants);

% --- STAGE 5: Neutral Density and Collision Profiles ---
neutral_profile = compute_neutral_density_profile(r_fine, constants);

% --- Verification Step: Plotting All Key Profiles ---
fprintf('Verification: Plotting all key analysis and model comparison results...\n');
% The normalised radius is created here for plotting purposes ONLY.
r_norm = r_fine / constants.machine.a;

% % --- Figure 1: Velocity Profile (Matches Thesis Fig. 3.5) ---
% figure('Name', 'Velocity Profile Analysis');
% hold on; grid on; box on;
% errorbar(exp_data.r_Vphi_exp / constants.machine.a, exp_data.Vphi_exp, exp_data.Vphi_err_exp, 'ko', 'MarkerFaceColor', 'k', 'DisplayName', 'Experimental Data');
% plot(r_norm, velocity_results.poly_fit_avg / 1000, 'g--', 'LineWidth', 2.5, 'DisplayName', '5th-order polynomial fit');
% fill([r_norm; flipud(r_norm)], [velocity_results.poly_fit_ci_lower / 1000; flipud(velocity_results.poly_fit_ci_upper / 1000)], 'g', 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', '95% Confidence interval');
% plot(r_norm, velocity_results.bessel_fit / 1000, 'r-', 'LineWidth', 2, 'DisplayName', 'Fourier-Bessel series fit');
% xlabel('Normalised Radius (r/a)'); ylabel('Toroidal Velocity, V_{\phi} (km/s)');
% title('Velocity Profile Analysis'); xlim([0, 1]); legend('show', 'Location', 'best');
% hold off;
% 
% % --- Figure 2: Ion Temperature Profile (Matches Thesis Fig. 3.6) ---
% figure('Name', 'Ion Temperature Profile Analysis');
% hold on; grid on; box on;
% errorbar(exp_data.r_Ti_exp / constants.machine.a, exp_data.Ti_reconstructed_exp, exp_data.Ti_reconstructed_err_exp, 'ko', 'MarkerFaceColor', 'k', 'DisplayName', 'Experimental Data');
% plot(r_norm, temperature_results.profile_avg, 'r-', 'LineWidth', 2.5, 'DisplayName', 'Canonical fit');
% fill([r_norm; flipud(r_norm)], [temperature_results.profile_ci_lower; flipud(temperature_results.profile_ci_upper)], 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', '95% Confidence Band');
% xlabel('Normalised Radius (r/a)'); ylabel('Ion Temperature, T_{i} (eV)');
% title('Ion Temperature Profile Analysis'); xlim([0, 1]); legend('show', 'Location', 'best');
% hold off;
% 
% % --- Figure 3: Thermal Velocity Profiles (Matches Thesis Fig. B.2) ---
% figure('Name', 'Thermal Velocity Profiles');
% hold on; grid on; box on;
% yyaxis left;
% plot(r_norm, derived_profiles.v_th_e, '-', 'LineWidth', 2);
% ylabel('Electron Thermal Velocity, v_{th,e} (m/s)');
% ylim([0, 2e7]);
% ax = gca; ax.YColor = [0 0.4470 0.7410];
% yyaxis right;
% plot(r_norm, derived_profiles.v_th_i, '-', 'LineWidth', 2);
% ylabel('Ion Thermal Velocity, v_{th,i} (m/s)');
% ylim([0, 3e5]);
% ax = gca; ax.YColor = [0.8500 0.3250 0.0980];
% xlabel('Normalised Radius (r/a)');
% title('Thermal Velocity Profiles'); xlim([0, 1]);
% hold off;
% 
% % --- Figure 4: Coulomb Logarithm Profiles (Matches Thesis Fig. B.3) ---
% figure('Name', 'Coulomb Logarithm Profiles');
% hold on; grid on; box on;
% plot(r_norm, derived_profiles.ln_lambda_ee, 'LineWidth', 2, 'DisplayName', '\lambda_{ee}');
% plot(r_norm, derived_profiles.ln_lambda_ei, 'LineWidth', 2, 'DisplayName', '\lambda_{ei}');
% plot(r_norm, derived_profiles.ln_lambda_ii, 'LineWidth', 2, 'DisplayName', '\lambda_{ii}');
% xlabel('Normalised Radius (r/a)'); ylabel('Coulomb Logarithm, \lambda = ln\Lambda');
% title('Coulomb Logarithm Profiles'); legend('show', 'Location', 'best'); xlim([0, 1]);
% hold off;
% 
% % --- Figure 5: Collision Metrics (Matches Thesis Fig. B.4) ---
% figure('Name', 'Collision Metrics');
% subplot(2, 2, 1); hold on; grid on; box on; plot(r_norm, derived_profiles.tau_e_calc, '--'); plot(r_norm, derived_profiles.tau_e_15, '-'); plot(r_norm, derived_profiles.tau_e_17, '-.'); xlabel('r/a'); ylabel('\tau_e [s]'); title('Electron Collision Time'); legend('\tau_e (\lambda_{ei,calc})', '\tau_e (\lambda = 15)', '\tau_e (\lambda = 17)'); xlim([0, 1]);
% subplot(2, 2, 2); hold on; grid on; box on; plot(r_norm, derived_profiles.tau_i_calc, '--'); plot(r_norm, derived_profiles.tau_i_15, '-'); plot(r_norm, derived_profiles.tau_i_17, '-.'); xlabel('r/a'); ylabel('\tau_i [s]'); title('Ion Collision Time'); legend('\tau_i (\lambda_{ii,calc})', '\tau_i (\lambda = 15)', '\tau_i (\lambda = 17)'); xlim([0, 1]);
% subplot(2, 2, 3); hold on; grid on; box on; plot(r_norm, derived_profiles.nu_ei / 1e5, '-'); xlabel('r/a'); ylabel('\nu_e [\times10^5 s^{-1}]'); title('Electron Collision Frequency'); xlim([0, 1]);
% subplot(2, 2, 4); hold on; grid on; box on; plot(r_norm, derived_profiles.nu_i / 1e3, '-'); xlabel('r/a'); ylabel('\nu_i [\times10^3 s^{-1}]'); title('Ion Collision Frequency'); xlim([0, 1]);
% 
% % --- Figure 6: Key Dimensionless Profiles (Matches Thesis Fig. B.5) ---
% figure('Name', 'Key Dimensionless Profiles');
% subplot(1, 2, 1); hold on; grid on; box on; plot(r_norm, derived_profiles.q_profile, 'b-', 'LineWidth', 2); plot(r_norm, derived_profiles.s_profile, 'r--', 'LineWidth', 2); xlabel('r/a'); ylabel('Value'); title('q and s Profiles'); legend('Safety Factor q', 'Magnetic Shear s'); xlim([0, 1]); ylim([0, 4]);
% subplot(1, 2, 2); hold on; grid on; box on; plot(r_norm, derived_profiles.nu_star_i, 'b-', 'LineWidth', 2); plot(r_norm, derived_profiles.nu_star_e, 'r--', 'LineWidth', 2); xlabel('r/a'); ylabel('Collisionality, \nu_*'); title('Collisionality Profiles'); legend('Ion Collisionality (\nu_*)', 'Electron Collisionality (\nu_*)'); xlim([0, 1]); ylim([0, 4]);

% --- Figure A: Helander Model vs. Experiment (Matches Thesis Fig. 3.4) ---
figure('Name', 'Helander Model Comparison');
hold on; grid on; box on;
% Plot experimental data and polynomial fit from velocity analysis
h1 = errorbar(exp_data.r_Vphi_exp / constants.machine.a, ...
    -exp_data.Vphi_exp, exp_data.Vphi_err_exp, ...
    'ko', 'MarkerFaceColor', 'k', 'DisplayName', 'Experimental data');
h2 = fill([r_norm; flipud(r_norm)], [-velocity_results.poly_fit_ci_lower / 1000; ...
    -flipud(velocity_results.poly_fit_ci_upper / 1000)], ...
    'g', 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', ...
    '95% Confidence interval');
h3 = plot(r_norm, -velocity_results.poly_fit_avg / 1000, 'g--', 'LineWidth', 2.5, 'DisplayName', '5th-order polynomial fit');
% Plot the Helander model prediction
h4 = plot(r_norm, -theoretical_models.Vphi_Helander / 1000, 'r-', 'LineWidth', 2.5, 'DisplayName', 'Helander model');
xlabel('Normalised Radius (r/a)');
ylabel('Toroidal Velocity, V_{\phi} (km/s)');
title('Comparison with Helander Model');
legend([h1, h3, h2, h4], 'Location', 'best');
xlim([0, 1]); ylim([-10, 30]);
hold off;

% --- Figure B: Solomon Model Profiles (Matches Thesis Fig. 3.9 & 3.10) ---
figure('Name', 'Solomon Model Profiles');
r_norm = r_fine / constants.machine.a;

% Subplot 1: Momentum Diffusivity
subplot(1, 2, 1);
hold on; grid on; box on;

% Isolate the data structure for the R/Ln ~ 6.1 case to ensure consistency
chi_calc = theoretical_models.chi_phi_solo_Rln_calc;
label_calc_val = theoretical_models.R_over_Ln_mid_value;
label_calc_profile = sprintf('\\chi_{\\phi}^{(Solo)} (R/L_n \\approx %.1f)', label_calc_val);
label_calc_fill = sprintf('Uncertainty (R/L_n \\approx %.1f)', label_calc_val);

% Plot the uncertainty band for the R/Ln ~ 6.1 case using its own data
fill([r_norm; flipud(r_norm)], [chi_calc.lower_bound; flipud(chi_calc.upper_bound)], ...
    [0.8500 0.3250 0.0980], 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', label_calc_fill);

% Plot the main profile line for the R/Ln ~ 6.1 case
plot(r_norm, chi_calc.profile, 'r-', 'LineWidth', 2.5, 'DisplayName', label_calc_profile);

% Plot profile for R/Ln = 2
chi_2 = theoretical_models.chi_phi_solo_Rln2;
h_line2 = plot(r_norm, chi_2.profile, 'b--', 'LineWidth', 2, 'DisplayName', '\chi_{\phi}^{(Solo)} (R/L_n = 2)');

xlabel('Normalised Radius (r/a)');
ylabel('Momentum Diffusivity, \chi_{\phi} [m^2/s]');
title('Solomon Momentum Diffusivity');
xlim([0, 1]); ylim([0, 25]);
legend('show', 'Location', 'northwest');

% Subplot 2: Pinch Velocity
subplot(1, 2, 2);
hold on; grid on; box on;
v_pinch = theoretical_models.Vpinch_solo;
fill([r_norm; flipud(r_norm)], [v_pinch.lower_bound; flipud(v_pinch.upper_bound)], [0 0.4470 0.7410], 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', 'Uncertainty');
plot(r_norm, v_pinch.profile, 'b-', 'LineWidth', 2, 'DisplayName', 'V_{pinch}^{(Solo)}');
xlabel('Normalised Radius (r/a)');
ylabel('Pinch Velocity, V_{pinch} [m/s]');
title('Solomon Pinch Velocity');
xlim([0, 1]); ylim([-60, 0]);
legend('show', 'Location', 'southwest');

% --- Figure C: Neutral Density Profile (Matches Thesis Fig. 3.7) ---
figure('Name', 'Neutral Density Profile');
hold on; grid on; box on;
plot(r_norm, neutral_profile.n_H0 / 1e16, 'LineWidth', 3);
xlabel('Normalised Radius (r/a)');
ylabel('Neutral Density, n_{H0} (\times10^{16} m^{-3})');
title('Neutral Particle Density Profile');
xlim([0, 1]); ylim([0, 6]);
hold off;

fprintf('Script finished.\n');