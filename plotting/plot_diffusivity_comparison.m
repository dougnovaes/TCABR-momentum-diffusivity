function plot_diffusivity_comparison(effective_diffusivity_thesis, scaling_law_results, constants, r_fine)
%PLOT_DIFFUSIVITY_COMPARISON Plots the comparison of effective diffusivity profiles.
%   This function generates the key comparison plot, confronting the effective
%   diffusivity calculated via the thesis method with the various profiles
%   derived from theoretical scaling laws (Solomon, Hahm, Gürcan, Peeters).

fprintf('Generating final diffusivity comparison plot...\n');

r_norm = r_fine / constants.machine.a;

figure('Name', 'Comparison of Effective Diffusivity Models');
hold on; grid on; box on;

% --- 1. Plot the Effective Diffusivity from the Thesis Method ---
eff_thesis = effective_diffusivity_thesis;
h_fill = fill([r_norm; flipud(r_norm)], [eff_thesis.ci_lower; flipud(eff_thesis.ci_upper)], ...
    'c', 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'DisplayName', '95% CI (Thesis Method)');
h_main = plot(r_norm, eff_thesis.profile_avg, 'b-', 'LineWidth', 3, 'DisplayName', '\chi_{\phi, thesis}^{eff}');


% --- 2. Plot the Effective Diffusivity Profiles from Scaling Laws ---
h_solo = plot(r_norm, scaling_law_results.chi_eff_Solo, '--', 'LineWidth', 2, 'DisplayName', '\chi_{\phi, Solomon}^{eff}');
h_hahm = plot(r_norm, scaling_law_results.chi_eff_Hahm, ':', 'LineWidth', 2, 'DisplayName', '\chi_{\phi, Hahm}^{eff}');
h_gurc = plot(r_norm, scaling_law_results.chi_eff_Gurcan, '-.', 'LineWidth', 2, 'DisplayName', '\chi_{\phi, Gürcan}^{eff}');
h_peet = plot(r_norm, scaling_law_results.chi_eff_Peeters_Rln_calc, 'd-', 'MarkerSize', 3, 'DisplayName', '\chi_{\phi, Peeters}^{eff}');

% --- 3. Final Plot Formatting ---
xlabel('Normalised Radius (r/a)');
ylabel('Effective Diffusivity, \chi_{\phi}^{eff} [m^2/s]');
title('Comparison: Thesis Method vs. Scaling Law Models');
xlim([0, 1]); ylim([0, 25]);
legend('show', 'Location', 'northwest');
hold off;

fprintf('Final comparison plot generated.\n\n');

end