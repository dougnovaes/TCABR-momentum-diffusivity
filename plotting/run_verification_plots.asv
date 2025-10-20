function run_verification_plots(constants, exp_data, velocity_results, ...
    temperature_results, derived_profiles, theoretical_models, effective_diffusivity_thesis, ~, ~)
%RUN_VERIFICATION_PLOTS Generates all diagnostic and comparison plots.
%   This function centralises all plotting routines to visualise the results
%   from each stage of the main analysis script.

fprintf('Generating all verification plots...\n');

r_fine = exp_data.r_profiles;
r_norm = r_fine / constants.machine.a; % Normalised radius for plotting ONLY.

% --- Figure 1: Velocity Profile (Matches Thesis Fig. 3.5) ---
figure('Name', 'Velocity Profile Analysis');
hold on; grid on; box on;
errorbar(exp_data.r_Vphi_exp / constants.machine.a, -exp_data.Vphi_exp, exp_data.Vphi_err_exp, 'ko', 'MarkerFaceColor', 'k', 'DisplayName', 'Experimental Data');
plot(r_norm, -velocity_results.poly_fit_avg / 1000, 'g--', 'LineWidth', 2.5, 'DisplayName', '5th-order polynomial fit');
fill([r_norm; flipud(r_norm)], [-velocity_results.poly_fit_ci_lower / 1000; flipud(-velocity_results.poly_fit_ci_upper / 1000)], 'g', 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', '95% Confidence interval');
plot(r_norm, -velocity_results.bessel_fit / 1000, 'r-', 'LineWidth', 2, 'DisplayName', 'Fourier-Bessel series fit');
xlabel('Normalised Radius (r/a)'); ylabel('Toroidal Velocity, V_{\phi} (km/s)');
title('Velocity Profile Analysis'); xlim([0, 1]); legend('show', 'Location', 'best');
hold off;

% --- Figure 2: Ion Temperature Profile (Matches Thesis Fig. 3.6) ---
figure('Name', 'Ion Temperature Profile Analysis');
hold on; grid on; box on;
errorbar(exp_data.r_Ti_exp / constants.machine.a, exp_data.Ti_reconstructed_exp, exp_data.Ti_reconstructed_err_exp, 'ko', 'MarkerFaceColor', 'k', 'DisplayName', 'Reconstructed Data');
plot(r_norm, temperature_results.profile_avg, 'r-', 'LineWidth', 2.5, 'DisplayName', 'Canonical fit');
fill([r_norm; flipud(r_norm)], [temperature_results.profile_ci_lower; flipud(temperature_results.profile_ci_upper)], 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', '95% Confidence Band');
xlabel('Normalised Radius (r/a)'); ylabel('Ion Temperature, T_{i} (eV)');
title('Ion Temperature Profile Analysis'); xlim([0, 1]); legend('show', 'Location', 'best');
hold off;

% --- Figure 3: Key Dimensionless Profiles (Matches Thesis Fig. B.5) ---
figure('Name', 'Key Dimensionless Profiles');
subplot(1, 2, 1); hold on; grid on; box on; plot(r_norm, derived_profiles.q_profile, 'b-', 'LineWidth', 2); plot(r_norm, derived_profiles.s_profile, 'r--', 'LineWidth', 2); xlabel('r/a'); ylabel('Value'); title('q and s Profiles'); legend('Safety Factor q', 'Magnetic Shear s'); xlim([0, 1]); ylim([0, 4]);
subplot(1, 2, 2); hold on; grid on; box on; plot(r_norm, derived_profiles.nu_star_i, 'b-', 'LineWidth', 2); plot(r_norm, derived_profiles.nu_star_e, 'r--', 'LineWidth', 2); xlabel('r/a'); ylabel('Collisionality, \nu_*'); title('Collisionality Profiles'); legend('Ion Collisionality (\nu_*)', 'Electron Collisionality (\nu_*)'); xlim([0, 1]); ylim([0, 4]);

% --- Figure 4: Helander Model vs. Experiment (Matches Thesis Fig. 3.4) ---
figure('Name', 'Helander Model Comparison');
hold on; grid on; box on;
% NOTE: Experimental V_phi is counter-current (negative). The Helander model
% as derived gives a co-current (positive) result. We negate the model's
% result for direct visual comparison on the same axes.
h1 = errorbar(exp_data.r_Vphi_exp / constants.machine.a, -exp_data.Vphi_exp, exp_data.Vphi_err_exp, 'ko', 'MarkerFaceColor', 'k', 'DisplayName', 'Experimental data');
h2 = fill([r_norm; flipud(r_norm)], [-velocity_results.poly_fit_ci_lower / 1000; flipud(-velocity_results.poly_fit_ci_upper / 1000)], 'g', 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', '95% Confidence interval');
h3 = plot(r_norm, -velocity_results.poly_fit_avg / 1000, 'g--', 'LineWidth', 2.5, 'DisplayName', '5th-order polynomial fit');
h4 = plot(r_norm, -theoretical_models.Vphi_Helander / 1000, 'r-', 'LineWidth', 2.5, 'DisplayName', 'Helander model');
xlabel('Normalised Radius (r/a)'); ylabel('Toroidal Velocity, V_{\phi} (km/s)');
title('Comparison with Helander Model');
legend([h1, h3, h2, h4], 'Location', 'best');
xlim([0, 1]); ylim([-10, 30]);
hold off;

% --- Figure 5: Effective Momentum Diffusivity (Matches Thesis Fig. 3.8) ---
figure('Name', 'Effective Momentum Diffusivity');
hold on; grid on; box on;
fill([r_norm; flipud(r_norm)], [effective_diffusivity_thesis.ci_lower; flipud(effective_diffusivity_thesis.ci_upper)], 'c', 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'DisplayName', '95% Confidence Band');
plot(r_norm, effective_diffusivity_thesis.profile_avg, 'b-', 'LineWidth', 2.5, 'DisplayName', 'Mean Effective Diffusivity');
xlabel('Normalised Radius (r/a)');
ylabel('Effective Diffusivity, \chi_{\phi}^{eff} [m^2/s]');
title('Effective Momentum Diffusivity (Thesis Method)');
xlim([0, 1]); ylim([0, 10]);
legend('show', 'Location', 'best');
hold off;

fprintf('Verification plots generated.\n\n');
end