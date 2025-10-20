function plot_stage1_justification_and_fits(constants, exp_data, velocity_results, temperature_results, theoretical_models, r_fine)
%PLOT_STAGE1_JUSTIFICATION_AND_FITS Generates plots for the primary fits and theoretical justification.
%   This function creates three key figures that replicate the narrative of the thesis:
%   1. The canonical fit of the ion temperature profile.
%   2. The comparison of the Helander model against experimental velocity data, which
%      justifies the theoretical approach of the thesis.
%   3. The full analysis of the velocity profile, including the Fourier-Bessel fit.

fprintf('Plotting Stage 1: Justification and Primary Fits...\n');

r_norm = r_fine / constants.machine.a;

% --- Figure 1: Ion Temperature Profile Analysis (Matches Thesis Fig. 3.6) ---
fig1 = figure('Name', 'Ion Temperature Profile Analysis');
ax1 = gca;
hold(ax1, 'on');

errorbar(ax1, exp_data.r_Ti_exp / constants.machine.a, exp_data.Ti_reconstructed_exp, ...
    exp_data.Ti_reconstructed_err_exp, 'ko', 'MarkerFaceColor', 'k', ...
    'MarkerSize', 8, 'LineWidth', 1.5, 'DisplayName', 'Experimental Data');
plot(ax1, r_norm, temperature_results.profile_avg, 'r-', 'LineWidth', 4, 'DisplayName', 'Canonical fit');
fill(ax1, [r_norm; flipud(r_norm)], [temperature_results.profile_ci_lower; flipud(temperature_results.profile_ci_upper)], ...
    'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', '95\% Confidence Band');
    
xlabel(ax1, 'Normalised Radius ($r/a$)');
ylabel(ax1, 'Ion Temperature, $T_{i}$ (eV)');
title(ax1, 'Ion Temperature Profile Analysis');
xlim(ax1, [0, 1]);
legend(ax1, 'show', 'Location', 'northeast');
set_publication_style(ax1); % Apply standard style
hold(ax1, 'off');


% --- Figure 2: Helander Model Justification (Matches Thesis Fig. 3.4) ---
fig2 = figure('Name', 'Helander Model Justification');
ax2 = gca;
hold(ax2, 'on');

% NOTE: Signs are flipped for visual comparison, as in the thesis figure.
h1 = errorbar(ax2, exp_data.r_Vphi_exp / constants.machine.a, -exp_data.Vphi_exp, ...
    exp_data.Vphi_err_exp, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8, ...
    'LineWidth', 1.5, 'DisplayName', 'Experimental data');
h2 = fill(ax2, [r_norm; flipud(r_norm)], [-velocity_results.poly_fit_ci_upper / 1000; -flipud(velocity_results.poly_fit_ci_lower / 1000)], ...
    'g', 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', '95\% Confidence interval');
h3 = plot(ax2, r_norm, -velocity_results.poly_fit_avg / 1000, 'g--', 'LineWidth', 5, 'DisplayName', '5th-order polynomial fit');
h4 = plot(ax2, r_norm, -theoretical_models.Vphi_Helander / 1000, 'r-', 'LineWidth', 4, 'DisplayName', 'Helander model');

yline(ax2, 0, '--', 'HandleVisibility', 'off'); % Added y=0 line

xlabel(ax2, 'Normalised Radius ($r/a$)');
ylabel(ax2, 'Toroidal Velocity, $V_{\phi}$ (km/s)');
title(ax2, 'Comparison with Helander Model');
legend(ax2, [h1, h3, h2, h4], 'Location', 'best');
xlim(ax2, [0, 1]); ylim(ax2, [-10, 30]);
set_publication_style(ax2); % Apply standard style
hold(ax2, 'off');


% --- Figure 3: Velocity Profile Analysis (Matches Thesis Fig. 3.5) ---
fig3 = figure('Name', 'Velocity Profile Analysis');
ax3 = gca;
hold(ax3, 'on');

% NOTE: Original (negative) velocity data is plotted here.
errorbar(ax3, exp_data.r_Vphi_exp / constants.machine.a, -exp_data.Vphi_exp, ...
    exp_data.Vphi_err_exp, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8, ...
    'LineWidth', 1.5, 'DisplayName', 'Experimental Data');
plot(ax3, r_norm, -velocity_results.poly_fit_avg / 1000, 'g--', 'LineWidth', 6, 'DisplayName', '5th-order polynomial fit');
fill(ax3, [r_norm; flipud(r_norm)], [-velocity_results.poly_fit_ci_lower / 1000; flipud(-velocity_results.poly_fit_ci_upper / 1000)], ...
    'g', 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', '95\% Confidence interval');
plot(ax3, r_norm, -velocity_results.bessel_fit / 1000, 'r-', 'LineWidth', 4, 'DisplayName', 'Fourier-Bessel series fit');

yline(ax3, 0, '--', 'HandleVisibility', 'off'); % Added y=0 line

xlabel(ax3, 'Normalised Radius ($r/a$)');
ylabel(ax3, 'Toroidal Velocity, $V_{\phi}$ (km/s)');
title(ax3, 'Velocity Profile Fourier-Bessel Analysis');
xlim(ax3, [0, 1]); ylim(ax3, [-10, 30]);
legend(ax3, 'show', 'Location', 'best');
set_publication_style(ax3); % Apply standard style
hold(ax3, 'off');

end