function plot_final_comparison_variants(diffusivity_exp, variants, constants, save_dir)
%PLOT_FINAL_COMPARISON_VARIANTS Generates a comprehensive set of experiment vs. theory comparison plots.
%   This function serves as the primary visualization tool for the final
%   analysis. It dynamically discovers all theoretical pinch models (e.g.,
%   Solomon, Hahm, Peeters) available in the 'variants' structure and
%   generates a separate, publication-quality figure for each one.
%
%   For each model, the plot confronts the experimentally-derived effective
%   diffusivity (with its 95% CI) against the four key theoretical variants,
%   each also including their Monte Carlo-derived 95% CIs:
%     1. Global R0, Mid-radius gradient
%     2. Global R0, Full radial gradient profile
%     3. Local R(r), Mid-radius gradient
%     4. Local R(r), Full radial gradient profile

    arguments
        diffusivity_exp (1,1) struct
        variants (1,1) struct
        constants (1,1) struct
        save_dir (1,1) string
    end
    
    fprintf('Generating final comparison plots with all theoretical variants...\n');
    
    % --- 1. Setup and Style Definitions ---
    r_norm = variants.meta.r_norm;
    colors = get(groot,'DefaultAxesColorOrder');
    
    % Define consistent styles for the four variants
    variant_styles = {
        {'DisplayName', 'Theory ($R_0$, grad @ mid)',    'LineStyle', '--', 'Color', colors(2,:), 'LineWidth', 2.5},
        {'DisplayName', 'Theory ($R_0$, grad profile)',  'LineStyle', '-.', 'Color', colors(3,:), 'LineWidth', 2.0},
        {'DisplayName', 'Theory ($R(r)$, grad @ mid)',   'LineStyle', ':',  'Color', colors(4,:), 'LineWidth', 2.5},
        {'DisplayName', 'Theory ($R(r)$, grad profile)', 'LineStyle', '--', 'Color', colors(5,:), 'LineWidth', 2.0}
    };
    variant_keys = {'R0_mid', 'R0_profile', 'Rlocal_mid', 'Rlocal_profile'};

    % Get a list of all pinch models calculated in the variants struct
    pinch_models = fieldnames(variants);
    pinch_models = pinch_models(~strcmp(pinch_models, 'meta')); % Exclude metadata field

    % --- 2. Loop Through Each Model and Generate a Figure ---
    for i = 1:numel(pinch_models)
        model_name = pinch_models{i};
        model_data = variants.(model_name);
        
        fig = figure('Name', ['Comparison vs. ', model_name], 'WindowStyle', 'docked');
        ax = gca;
        hold(ax, 'on');

        % --- 2a. Plot Experimental Data as the Reference ---
        fill(ax, [r_norm; flipud(r_norm)], ...
            [diffusivity_exp.ci_lower; flipud(diffusivity_exp.ci_upper)], ...
            colors(1,:), 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'DisplayName', '95% CI (Exp.)');
            
        plot(ax, r_norm, diffusivity_exp.profile_avg, '-', 'LineWidth', 3.5, ...
            'Color', colors(1,:), 'DisplayName', '$\chi_{\phi, \mathrm{eff}}^{(\mathrm{exp})}$');

        % --- 2b. Plot the Pure Diffusivity Baseline ---
        plot(ax, r_norm, variants.meta.chi_phi_solo_avg, ':', 'LineWidth', 2, 'Color', 'k', ...
            'DisplayName', '$\chi_{\phi}^{(\mathrm{Solo})}$ (pure)');

        % --- 2c. Loop Through and Plot the Four Variants for the Current Model ---
        for j = 1:numel(variant_keys)
            key = variant_keys{j};
            if isfield(model_data, key)
                var_data = model_data.(key);
                style = variant_styles{j};
                
                % Plot confidence interval fill
                fill(ax, [r_norm; flipud(r_norm)], [var_data.ci_lower; flipud(var_data.ci_upper)], ...
                    style{6}, 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'HandleVisibility', 'off');
                
                % Plot mean profile line
                plot(ax, r_norm, var_data.profile_avg, style{:});
            end
        end

        % --- 2d. Formatting and Annotations ---
        xlabel('Normalised Radius ($r/a$)');
        ylabel('Effective Diffusivity, $\chi_{\phi, \mathrm{eff}}$ [m$^2$/s]');
        title(sprintf('Effective Momentum Diffusivity: Experiment vs. %s Model', ...
            strrep(model_name, '_', '\_')));
        xlim([0, 1]);
        ylim_upper = max(diffusivity_exp.ci_upper(r_norm > 0.1 & r_norm < 0.9));
        ylim([0, ylim_upper * 1.8]); % Dynamic y-limit
        legend('show', 'Location', 'northwest');
        set_publication_style(ax);

        % Add annotation with key parameters
        anno_text = sprintf('$R/L_{n, \\mathrm{mid}} = %.2f$\n$R/L_{V\\phi, \\mathrm{mid}} = %.2f$', ...
            variants.meta.R_over_Ln_mid, variants.meta.R_over_LVphi_mid);
        annotation('textbox', [0.65, 0.75, 0.25, 0.15], 'String', anno_text, ...
            'Interpreter', 'latex', 'FontSize', 18, 'FitBoxToText', 'on', ...
            'BackgroundColor', 'white', 'EdgeColor', 'k');
        
        hold(ax, 'off');
        
        % --- 2e. Save the Figure ---
        fig_filename = fullfile(save_dir, ['final_comparison_', model_name, '.png']);
        saveas(fig, fig_filename);
        savefig(fig, strrep(fig_filename, '.png', '.fig'));
        fprintf('...plot saved to %s\n', fig_filename);
    end
end