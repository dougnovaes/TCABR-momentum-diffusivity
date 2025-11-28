function plot_advanced_fit_results(results_adv, save_dir)
%PLOT_ADVANCED_FIT_RESULTS Visualizes results of the advanced inverse analysis.

    arguments
        results_adv (1,1) struct
        save_dir (1,1) string
    end
    
    fprintf('Generating advanced fit analysis plots...\n');
    
    r_norm = results_adv.r_norm;
    target = results_adv.target;
    fit = results_adv.chi_fit;
    nH0 = results_adv.nH0_inferred;
    mask = results_adv.mask;
    params = results_adv.p_opt;
    
    adv_dir = fullfile(save_dir, 'advanced_studies');
    if ~exist(adv_dir, 'dir'), mkdir(adv_dir); end

    % Plot 1: Diffusivity Fit
    f1 = figure('Name', 'Advanced Fit: Diffusivity', 'WindowStyle', 'docked');
    ax1 = gca; hold on; box on; grid on;
    
    % Shaded excluded regions
    y_max = max(max(target(isfinite(target))), max(fit));
    area(r_norm, ~mask * y_max*1.2, 'FaceColor', [0.9 0.9 0.9], 'EdgeColor', 'none', 'DisplayName', 'Excluded Region');
    
    plot(ax1, r_norm, target, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Benchmark: Solomon (Local)');
    plot(ax1, r_norm, fit, 'r-', 'LineWidth', 2.5, 'DisplayName', 'Optimised Fit (Voigt)');
    
    xlabel(ax1, 'r/a'); ylabel(ax1, '\chi_{\phi} [m^2/s]');
    title(ax1, 'Inverse Analysis: Fitting Effective Diffusivity');
    legend(ax1, 'Location', 'northwest');
    ylim(ax1, [0, 65]); xlim(ax1, [0, 1]);
    
    txt = sprintf('Offset: %.2f m^2/s\nCenter: %.2f\nWidth: %.3f', params(5), params(2), params(3));
    text(ax1, 0.05, 0.5, txt, 'Units', 'normalized', 'BackgroundColor', 'w', 'EdgeColor', 'k');
    
    saveas(f1, fullfile(adv_dir, 'adv_fit_diffusivity.png'));
    
    % Plot 2: Inferred Neutrals
    f2 = figure('Name', 'Inferred Neutral Density', 'WindowStyle', 'docked');
    ax2 = gca; hold on; box on; grid on;
    plot(ax2, r_norm, nH0, 'b-', 'LineWidth', 2.5, 'DisplayName', 'Inferred n_{H0}');
    xlabel(ax2, 'r/a'); ylabel(ax2, 'n_{H0} [m^{-3}]');
    title(ax2, 'Inferred Neutral Density Profile');
    legend(ax2, 'Location', 'best');
    xlim(ax2, [0, 1]);
    
    saveas(f2, fullfile(adv_dir, 'adv_fit_neutrals.png'));
    fprintf('...advanced plots saved to %s\n', adv_dir);
end