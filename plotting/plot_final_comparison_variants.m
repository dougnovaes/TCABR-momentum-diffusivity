function plot_final_comparison_variants(diffusivity_exp, variants, constants, save_dir)
%PLOT_FINAL_COMPARISON_VARIANTS Generates a comprehensive set of experiment vs. theory comparison plots.
%   This function creates a separate, publication-quality figure for each
%   theoretical model found in the 'variants' structure.
%
%   Each figure includes:
%   1. A "Theoretical Annotation Box" displaying the relevant physics formulas.
%   2. An "intelligent" legend that incorporates the numerical values of the
%      key parameters (`R/L_n,mid` and `R/L_Vphi,mid`) directly into the legend text.
%   3. A detailed comparison of the experimental data against the four main
%      theoretical variants for that model, all with 95% confidence intervals.

    arguments
        diffusivity_exp (1,1) struct
        variants (1,1) struct
        constants (1,1) struct
        save_dir (1,1) string
    end
    
    fprintf('Generating final comparison plots with formula annotations...\n');
    
    % --- 1. Setup, Styles, and Formula Definitions ---
    r_norm = variants.meta.r_norm;
    colors = get(groot,'DefaultAxesColorOrder');
    meta = variants.meta;

    % Create a map of LaTeX formulas for each pinch model
    formula_map = containers.Map('KeyType', 'char', 'ValueType', 'char');
    formula_map('Solomon') = '$V_{\mathrm{pinch}} = C \cdot \nu_e^*$';
    formula_map('Hahm') = '$V_{\mathrm{pinch}} = -2 \chi_{\phi} / R_0$';
    formula_map('Gurcan') = '$V_{\mathrm{pinch}} = -(2\chi_{\phi}/R) \cdot (F+r/R_0)$';
    formula_map('Peeters_Rln2') = '$V_{\mathrm{pinch}} = (\chi_{\phi}/R) \cdot (-4 - R/L_n)$, with $R/L_n=2$';
    formula_map('Peeters_Rln_calc') = '$V_{\mathrm{pinch}} = (\chi_{\phi}/R) \cdot (-4 - R/L_n)$';
    
    % General formula for all plots
    general_formula = '$$\chi_{\phi, \mathrm{eff}} = \chi_{\phi}^{(\mathrm{Solo})} \left(1 + \frac{R \cdot V_{\mathrm{pinch}}}{\chi_{\phi}^{(\mathrm{Solo})}} \frac{1}{R/L_{V_{\phi}}}\right)$$';
    
    % Define styles and legend text formats
    variant_styles = {
        {'DisplayName', 'Theory ($R_0$, grad @ mid = %.2f)',    'LineStyle', '--', 'Color', colors(2,:), 'LineWidth', 2.5},
        {'DisplayName', 'Theory ($R_0$, grad profile)',           'LineStyle', '-.', 'Color', colors(3,:), 'LineWidth', 2.0},
        {'DisplayName', 'Theory ($R(r)$, grad @ mid = %.2f)',   'LineStyle', ':',  'Color', colors(4,:), 'LineWidth', 2.5},
        {'DisplayName', 'Theory ($R(r)$, grad profile)',          'LineStyle', '--', 'Color', colors(5,:), 'LineWidth', 2.0}
    };
    variant_keys = {'R0_mid', 'R0_profile', 'Rlocal_mid', 'Rlocal_profile'};

    pinch_models = fieldnames(variants);
    pinch_models = pinch_models(~strcmp(pinch_models, 'meta'));

    % --- 2. Loop Through Each Model and Generate a Figure ---
    for i = 1:numel(pinch_models)
        model_name = pinch_models{i};
        model_data = variants.(model_name);
        
        fig = figure('Name', ['Comparison vs. ', model_name], 'WindowStyle', 'docked', 'Position', [100, 100, 1100, 800]);
        ax = gca;
        hold(ax, 'on');

        % --- 2a. Plot Experimental Data ---
        fill(ax, [r_norm; flipud(r_norm)], [diffusivity_exp.ci_lower; flipud(diffusivity_exp.ci_upper)], ...
            colors(1,:), 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'DisplayName', '95% CI (Exp.)');
        plot(ax, r_norm, diffusivity_exp.profile_avg, '-', 'LineWidth', 3.5, ...
            'Color', colors(1,:), 'DisplayName', '$\chi_{\phi, \mathrm{eff}}^{(\mathrm{exp})}$');

        % --- 2b. Plot Pure Diffusivity Baseline ---
        plot(ax, r_norm, meta.chi_phi_solo_avg, ':', 'LineWidth', 2, 'Color', 'k', ...
            'DisplayName', '$\chi_{\phi}^{(\mathrm{Solo})}$ (pure)');

        % --- 2c. Loop Through and Plot the Four Variants ---
        plot_handles = gobjects(numel(variant_keys), 1);
        for j = 1:numel(variant_keys)
            key = variant_keys{j};
            if isfield(model_data, key)
                var_data = model_data.(key);
                style = variant_styles{j};
                
                fill(ax, [r_norm; flipud(r_norm)], [var_data.ci_lower; flipud(var_data.ci_upper)], ...
                    style{6}, 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'HandleVisibility', 'off');
                
                % Format display name with value if it's a 'mid' variant
                if contains(key, 'mid')
                    display_name = sprintf(style{2}, meta.R_over_LVphi_mid);
                else
                    display_name = style{2};
                end
                
                plot_handles(j) = plot(ax, r_norm, var_data.profile_avg, style{1}, display_name, style{3:end});
            end
        end

        % --- 2d. Formatting, Legend and Annotations ---
        xlabel('Normalised Radius ($r/a$)');
        ylabel('Effective Diffusivity, $\chi_{\phi, \mathrm{eff}}$ [m$^2$/s]');
        title_str = sprintf('Effective Momentum Diffusivity: Experiment vs. %s Model', strrep(model_name, '_', '\_'));
        title(title_str);
        
        ylim_upper = max(diffusivity_exp.ci_upper(r_norm > 0.1 & r_norm < 0.9));
        xlim([0, 1]); ylim([0, ylim_upper * 1.8]);
        
        % Create the intelligent legend
        legend_title = sprintf('Theoretical Variants\n(Base $\\chi_\\phi$ uses $R/L_{n,mid}=%.2f$)', meta.R_over_Ln_mid);
        lgd = legend(ax, 'Location', 'northeast');
        lgd.Title.String = legend_title;
        lgd.Title.Interpreter = 'latex';
        
        set_publication_style(ax);
        
        % Create the theoretical formula annotation box
        model_formula = formula_map(model_name);
        anno_text = {general_formula, '', '\bf{Pinch Model:}', model_formula};
        annotation('textbox', [0.15, 0.6, 0.3, 0.3], 'String', anno_text, ...
            'Interpreter', 'latex', 'FontSize', 16, 'FitBoxToText', 'on', ...
            'BackgroundColor', [1, 1, 1, 0.8], 'EdgeColor', 'k');
        
        hold(ax, 'off');
        
        % --- 2e. Save the Figure ---
        fig_filename = fullfile(save_dir, ['final_comparison_', model_name, '.png']);
        saveas(fig, fig_filename);
        savefig(fig, strrep(fig_filename, '.png', '.fig'));
        fprintf('...plot saved to %s\n', fig_filename);
    end
end