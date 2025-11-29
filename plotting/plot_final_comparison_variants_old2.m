function plot_final_comparison_variants(diffusivity_exp, variants, constants, save_dir)
%PLOT_FINAL_COMPARISON_VARIANTS Generates a final set of publication-quality comparison plots.
%   This function creates a separate, highly polished figure for each theoretical model.
%
%   CHANGES:
%   - Plots ONLY the best theoretical variant: R=R0 (Global) with Local Vphi Gradients.
%   - Simplifies legends for clarity.
%   - Maintains robust LaTeX formula annotations and fixed Y-axis scaling.

    arguments
        diffusivity_exp (1,1) struct
        variants (1,1) struct
        constants (1,1) struct
        save_dir (1,1) string
    end
    
    fprintf('Generating final, publication-quality comparison plots...\n');
    
    % --- 1. Setup & Definitions ---
    r_norm = variants.meta.r_norm;
    colors = get(groot,'DefaultAxesColorOrder');
    meta = variants.meta;

    % --- LaTeX Formula Definitions (Preserving Layout) ---
    formula_map = containers.Map('KeyType', 'char', 'ValueType', 'any');
    
    % Base Diffusivity Formula
    % Note: Using \( \) for inline math in the annotation box context
    chi_solo_base = '\(\chi_{\varphi}^{(\mathrm{Solomon})} = (6.09 \pm 0.72)\nu_e^* + (0.157 \pm 0.072)R/L_n\)';
    rln_mid_text = sprintf(',\nwith \\(R/L_{n,mid}=%.2f\\)', meta.R_over_Ln_mid);
    chi_solo_formula = [chi_solo_base, rln_mid_text]; 

    % Pinch Velocity Formulas per Model
    formula_map('Solomon') = {chi_solo_formula, '\(V_{\mathrm{pinch}}^{(\mathrm{Solomon})} = (-24.2 \pm 3.5)\nu_e^*\)'};
    formula_map('Hahm')    = {chi_solo_formula, '\(V_{\mathrm{pinch}}^{(\mathrm{Hahm})} = -2 \chi_{\varphi} / R_0\)'};
    formula_map('Gurcan')  = {chi_solo_formula, '\(V_{\mathrm{pinch}}^{(\mathrm{G\ddot{u}rcan})} = -(2\chi_{\varphi}/R) \cdot (F+r/R_0)\)'};
    % Note: Peeters R/Ln=2 case is handled specifically in loop if needed, but formula is general
    formula_map('Peeters_Rln2') = {chi_solo_formula, '\(V_{\mathrm{pinch}}^{(\mathrm{Peeters})} = (\chi_{\varphi}/R) \cdot (-4 - R/L_n)\), with \(R/L_n=2\)'};
    formula_map('Peeters_Rln_calc') = {chi_solo_formula, '\(V_{\mathrm{pinch}}^{(\mathrm{Peeters})} = (\chi_{\varphi}/R) \cdot (-4 - R/L_n)\)'};
    
    % General Effective Diffusivity Framework
    general_formula = ['$$\chi_{\varphi}^{\mathrm{eff}} = \chi_{\varphi} \left(1 + \frac{R \cdot V_{\mathrm{pinch}}}{\chi_{\varphi}} ' ...
        '\frac{1}{R/L_{V_{\varphi}}}\right)$$'];

    % Identify models to plot
    pinch_models = fieldnames(variants);
    pinch_models = pinch_models(~strcmp(pinch_models, 'meta'));

    % --- 2. Loop Through Each Model and Generate Figure ---
    for i = 1:numel(pinch_models)
        model_name = pinch_models{i};
        model_data = variants.(model_name);
        
        fig = figure('Name', ['Comparison vs. ', model_name], 'WindowStyle', 'docked');
        ax = gca; hold(ax, 'on');

        % --- 2a. Plot Experimental Data (Blue) ---
        % Use experimental grid if available to match size
        if isfield(diffusivity_exp, 'r_fine')
            r_exp_norm = diffusivity_exp.r_fine / constants.machine.a;
        else
            r_exp_norm = r_norm;
        end
        
        fill(ax, [r_exp_norm; flipud(r_exp_norm)], [diffusivity_exp.ci_lower; flipud(diffusivity_exp.ci_upper)], ...
            colors(1,:), 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'DisplayName', '95\% CI (Exp.)');
        plot(ax, r_exp_norm, diffusivity_exp.profile_avg, '-', 'LineWidth', 3.5, ...
            'Color', colors(1,:), 'DisplayName', '$\chi_\varphi^{\mathrm{eff}}$ (Exp.)');

        % --- 2b. Plot Theoretical Baseline (Black Dotted) ---
        % This represents the turbulence level without pinch flux
        plot(ax, r_norm, meta.chi_phi_solo_avg, ':', 'LineWidth', 2.5, 'Color', 'k', ...
            'DisplayName', '$\chi_\varphi$ (Base/Solomon)');

        % --- 2c. Plot Theoretical Variant (Green Dashed) ---
        % We selected: R0_profile (Global R0, Local Vphi Gradient)
        target_key = 'R0_profile';
        
        if isfield(model_data, target_key)
            var_data = model_data.(target_key);
            
            % Theoretical Confidence Interval
            fill(ax, [r_norm; flipud(r_norm)], [var_data.ci_lower; flipud(var_data.ci_upper)], ...
                 colors(5,:), 'FaceAlpha', 0.10, 'EdgeColor', 'none', 'HandleVisibility', 'off');
            
            % Theoretical Mean Profile
            plot(ax, r_norm, var_data.profile_avg, '--', 'LineWidth', 3.0, ...
                'Color', colors(5,:), 'DisplayName', '$\chi_{\varphi}^{\mathrm{eff}}$ (Theory w/ Pinch)');
        else
            warning('Variant %s not found for model %s', target_key, model_name);
        end

        % --- 2d. Formatting ---
        xlabel('Normalised Radius ($r/a$)', 'Interpreter', 'latex');
        ylabel('Effective Diffusivity, $\chi_{\varphi}^{\mathrm{eff}}$ [m$^2$/s]', 'Interpreter', 'latex');
        
        % Clean Title
        clean_name = strrep(model_name, '_', '\_');
        if contains(model_name, 'Peeters')
             title_str = sprintf('Experiment vs. %s', clean_name);
        else
             title_str = sprintf('Experiment vs. %s Scaling', clean_name);
        end
        title(title_str, 'Interpreter', 'latex');
        
        % Axis Limits
        xlim([0, 1]); 
        ylim([0, 25]); % Fixed limit as requested
        
        % Legend (Simplified)
        lgd = legend(ax, 'Location', 'north');
        set(lgd, 'Interpreter', 'latex', 'Box', 'off', 'FontSize', 14); 
        
        set_publication_style(ax);
        
        % --- 2e. Annotation Box (Equations) ---
        if formula_map.isKey(model_name)
            model_formulas = formula_map(model_name);
            anno_text = { ...
                '\bf{Framework:}', general_formula, ...
                '\bf{Base Diffusivity:}', model_formulas{1}, ...
                '\bf{Pinch Model:}', model_formulas{2} ...
            };
            
            % Position box to avoid overlapping the curve peak if possible
            annotation('textbox', [0.20, 0.45, 0.4, 0.2], 'String', anno_text, ...
                'Interpreter', 'latex', 'FontSize', 14, 'FitBoxToText', 'on', ...
                'BackgroundColor', [1, 1, 1, 0.85], 'EdgeColor', 'k', 'VerticalAlignment', 'top');
        end
        
        hold(ax, 'off');
        
        % --- 2f. Save ---
        fig_filename = fullfile(save_dir, ['final_comparison_', model_name, '.png']);
        saveas(fig, fig_filename);
        savefig(fig, strrep(fig_filename, '.png', '.fig'));
        fprintf('...plot saved to %s\n', fig_filename);
    end
end