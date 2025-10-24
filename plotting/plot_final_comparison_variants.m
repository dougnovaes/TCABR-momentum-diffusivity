function plot_final_comparison_variants(diffusivity_exp, variants, constants, save_dir)
%PLOT_FINAL_COMPARISON_VARIANTS Generates a final set of publication-quality comparison plots.
%   This function creates a separate, highly polished figure for each
%   theoretical model. Each figure is richly annotated with theoretical
%   formulas and data-driven legends for maximum clarity and impact.

    arguments
        diffusivity_exp (1,1) struct
        variants (1,1) struct
        constants (1,1) struct
        save_dir (1,1) string
    end
    
    fprintf('Generating final, publication-quality comparison plots...\n');
    
    % --- 1. Setup, Styles, and Formula Definitions ---
    r_norm = variants.meta.r_norm;
    colors = get(groot,'DefaultAxesColorOrder');
    meta = variants.meta;

    % Create a map of LaTeX formulas for each pinch model.
    formula_map = containers.Map('KeyType', 'char', 'ValueType', 'any');
    chi_solo_formula = sprintf('$\\chi_{\\phi}^{(\\mathrm{Solo})} = (6.09 \\pm 0.72)\\nu_e^* + (0.157 \\pm 0.072)R/L_n$,\nwith $R/L_{n,mid}=%.2f$', meta.R_over_Ln_mid);
    formula_map('Solomon') = {chi_solo_formula, '$V_{\mathrm{pinch}}^{(\mathrm{Solo})} = (-24.2 \pm 3.5)\nu_e^*$';};
    formula_map('Hahm')    = {chi_solo_formula, '$V_{\mathrm{pinch}}^{(\mathrm{Hahm})} = -2 \chi_{\phi} / R_0$';};
    formula_map('Gurcan')  = {chi_solo_formula, '$V_{\mathrm{pinch}}^{(\mathrm{G\ddot{u}rc})} = -(2\chi_{\phi}/R) \cdot (F+r/R_0)$';};
    formula_map('Peeters_Rln2') = {chi_solo_formula, '$V_{\mathrm{pinch}}^{(\mathrm{Peet})} = (\chi_{\phi}/R) \cdot (-4 - R/L_n)$, with $R/L_n=2$';};
    formula_map('Peeters_Rln_calc') = {chi_solo_formula, '$V_{\mathrm{pinch}}^{(\mathrm{Peet})} = (\chi_{\phi}/R) \cdot (-4 - R/L_n)$';};
    
    general_formula = '$$\chi_{\phi, \mathrm{eff}} = \chi_{\phi} \left(1 + \frac{R \cdot V_{\mathrm{pinch}}}{\chi_{\phi}} \frac{1}{R/L_{V_{\phi}}}\right)$$';
    
    % Define styles and new, explicit legend text formats
    variant_styles = {
        {'DisplayName', '$R=R_0$, $R/L_{V_\phi}=%.2f$',    'LineStyle', '--', 'Color', colors(2,:), 'LineWidth', 2.5},
        {'DisplayName', '$R=R_0$, $R/L_{V_\phi}(r)$',      'LineStyle', '-.', 'Color', colors(3,:), 'LineWidth', 2.0},
        {'DisplayName', '$R=R(r)$, $R/L_{V_\phi}=%.2f$',   'LineStyle', ':',  'Color', colors(4,:), 'LineWidth', 2.5},
        {'DisplayName', '$R=R(r)$, $R/L_{V_\phi}(r)$',      'LineStyle', '--', 'Color', colors(5,:), 'LineWidth', 2.0}
    };
    variant_keys = {'R0_mid', 'R0_profile', 'Rlocal_mid', 'Rlocal_profile'};

    pinch_models = fieldnames(variants);
    pinch_models = pinch_models(~strcmp(pinch_models, 'meta'));

    % --- 2. Loop Through Each Model and Generate a Figure ---
    for i = 1:numel(pinch_models)
        model_name = pinch_models{i};
        model_data = variants.(model_name);
        
        fig = figure('Name', ['Comparison vs. ', model_name], 'WindowStyle', 'docked');
        ax = gca;
        hold(ax, 'on');

        % --- 2a. Plot Experimental Data ---
        fill(ax, [r_norm; flipud(r_norm)], [diffusivity_exp.ci_lower; flipud(diffusivity_exp.ci_upper)], ...
            colors(1,:), 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'DisplayName', '95\% CI (Exp.)');
        plot(ax, r_norm, diffusivity_exp.profile_avg, '-', 'LineWidth', 3.5, ...
            'Color', colors(1,:), 'DisplayName', '$\chi_{\phi, \mathrm{eff}}^{(\mathrm{exp})}$');

        % --- 2b. Plot Pure Diffusivity Baseline ---
        plot(ax, r_norm, meta.chi_phi_solo_avg, ':', 'LineWidth', 2, 'Color', 'k', ...
            'DisplayName', '$\chi_{\phi}^{(\mathrm{Solo})}$ (base)');

        % --- 2c. Loop Through and Plot the Four Variants ---
        for j = 1:numel(variant_keys)
            key = variant_keys{j};
            if isfield(model_data, key)
                var_data = model_data.(key);
                style = variant_styles{j};
                
                fill(ax, [r_norm; flipud(r_norm)], [var_data.ci_lower; flipud(var_data.ci_upper)], ...
                    style{6}, 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'HandleVisibility', 'off');
                
                if contains(key, 'mid')
                    display_name = sprintf(style{2}, meta.R_over_LVphi_mid);
                else
                    display_name = style{2};
                end
                
                plot(ax, r_norm, var_data.profile_avg, style{1}, display_name, style{3:end});
            end
        end

        % --- 2d. Formatting, Legend and Annotations ---
        xlabel('Normalised Radius ($r/a$)');
        ylabel('Effective Diffusivity, $\chi_{\phi, \mathrm{eff}}$ [m$^2$/s]');
        title_str = sprintf('Effective Momentum Diffusivity: Experiment vs. %s Model', strrep(model_name, '_', '\_'));
        title(title_str);
        
        ylim_upper = max(diffusivity_exp.ci_upper(r_norm > 0.1 & r_norm < 0.9));
        xlim([0, 1]); ylim([0, ylim_upper * 1.8]);
        
        % Create legend at the top-center
        lgd = legend(ax, 'Location', 'north');
        lgd.Title.String = 'Theoretical Variants'; % Simple title
        
        set_publication_style(ax);
        
        % Create the theoretical formula annotation box in the top-right
        model_formulas = formula_map(model_name);
        anno_text = { ...
            '\bf{Framework:}', general_formula, ...
            '\bf{Base Diffusivity:}', model_formulas{1}, ...
            '\bf{Pinch Model:}', model_formulas{2} ...
        };
        annotation('textbox', [0.55, 0.68, 0.4, 0.2], 'String', anno_text, ...
            'Interpreter', 'latex', 'FontSize', 13, 'FitBoxToText', 'on', ...
            'BackgroundColor', [1, 1, 1, 0.85], 'EdgeColor', 'k', 'VerticalAlignment', 'top');
        
        hold(ax, 'off');
        
        % --- 2e. Save the Figure ---
        fig_filename = fullfile(save_dir, ['final_comparison_', model_name, '.png']);
        saveas(fig, fig_filename);
        savefig(fig, strrep(fig_filename, '.png', '.fig'));
        fprintf('...plot saved to %s\n', fig_filename);
    end
end