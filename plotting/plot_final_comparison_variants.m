function plot_final_comparison_variants(diffusivity_exp, variants, save_dir)
%PLOT_FINAL_COMPARISON_VARIANTS Plots the definitive experimental vs. theoretical comparison.
%   This function visualizes the primary scientific results, overlaying the
%   experimentally derived effective diffusivity (with its 95% CI) against a
%   selection of the most relevant theoretical variants, which also include
%   their own Monte Carlo-derived 95% CIs.

    arguments
        diffusivity_exp (1,1) struct
        variants (1,1) struct
        save_dir (1,1) string
    end
    
    fprintf('Generating final comparison plot with theoretical variants...\n');
    
    fig = figure('Name', 'Final Comparison: Experimental vs. Theoretical Variants', 'Position', [100, 100, 1000, 750]);
    ax = gca;
    hold(ax, 'on');

    r_norm = variants.meta.r_norm;
    colors = get(groot,'DefaultAxesColorOrder');
    
    % --- 1. Plot Experimental Data (Thesis Method) ---
    fill(ax, [r_norm; flipud(r_norm)], ...
        [diffusivity_exp.ci_lower; flipud(diffusivity_exp.ci_upper)], ...
        colors(1,:), 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'DisplayName', '95% CI (Exp.)');
        
    plot(ax, r_norm, diffusivity_exp.profile_avg, '-', 'LineWidth', 3, ...
        'Color', colors(1,:), 'DisplayName', '$\chi_{\phi, \mathrm{eff}}^{(\mathrm{exp})}$');

    % --- 2. Plot a Selection of Theoretical Variants (e.g., for Solomon model) ---
    model_to_plot = 'Solomon';
    
    % Base Diffusivity
    plot(ax, r_norm, variants.meta.chi_phi_solo_avg, ':', 'LineWidth', 2, 'Color', 'k', ...
        'DisplayName', '$\chi_{\phi}^{(\mathrm{Solo})}$ (pure)');
        
    % Variant: R0, Mid-radius gradient
    var_data = variants.(model_to_plot).R0_mid;
    fill(ax, [r_norm; flipud(r_norm)], [var_data.ci_lower; flipud(var_data.ci_upper)], ...
        colors(2,:), 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'DisplayName', '95% CI (Theory)');
        
    plot(ax, r_norm, var_data.profile_avg, '--', 'LineWidth', 2.5, 'Color', colors(2,:), ...
        'DisplayName', sprintf('Theory (%s, $R_0$, grad @ mid)', model_to_plot));

    % Variant: R(r), Profile gradient
    var_data_2 = variants.(model_to_plot).Rlocal_profile;
    fill(ax, [r_norm; flipud(r_norm)], [var_data_2.ci_lower; flipud(var_data_2.ci_upper)], ...
        colors(3,:), 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off'); % This fill is just for shading, no legend entry
        
    plot(ax, r_norm, var_data_2.profile_avg, '-.', 'LineWidth', 2.5, 'Color', colors(3,:), ...
        'DisplayName', sprintf('Theory (%s, $R(r)$, grad profile)', model_to_plot));

    % --- 3. Formatting ---
    xlabel('Normalised Radius ($r/a$)');
    ylabel('Effective Diffusivity, $\chi_{\phi, \mathrm{eff}}$ [m$^2$/s]');
    title('Effective Momentum Diffusivity: Experiment vs. Theory');
    xlim([0, 1]);
    ylim_upper = max(diffusivity_exp.ci_upper(r_norm > 0.1));
    ylim([0, ylim_upper * 1.5]);
    legend('show', 'Location', 'northwest');
    set_publication_style(ax);

    % --- 4. Save Figure ---
    if ~exist(save_dir, 'dir'), mkdir(save_dir); end
    fig_filename = fullfile(save_dir, 'final_comparison_variants');
    saveas(fig, [fig_filename, '.png']);
    savefig(fig, [fig_filename, '.fig']);
    fprintf('...plot saved to %s.(png/fig)\n', fig_filename);
    hold(ax, 'off');
end