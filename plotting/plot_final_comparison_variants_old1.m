function plot_final_comparison_variants(diffusivity_exp, variants, constants, save_dir)
%PLOT_FINAL_COMPARISON_VARIANTS Generates a final set of publication-quality comparison plots.
%   This function creates a separate, highly polished figure for each
%   theoretical model. Each figure is richly annotated with theoretical
%   formulas and data-driven legends for maximum clarity and impact.
%
%   UPDATED: Now dynamically adjusts Y-axis limits to accommodate both
%   experimental data and theoretical curves, ignoring singularities.

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

    % Create a map of LaTeX formulas for each pinch model using robust delimiters.
    formula_map = containers.Map('KeyType', 'char', 'ValueType', 'any');
    
    % Construct chi_solo_formula carefully
    chi_solo_base = '\(\chi_{\varphi}^{(\mathrm{Solomon})} = (6.09 \pm 0.72)\nu_e^* + (0.157 \pm 0.072)R/L_n\)';
    rln_mid_text = sprintf(',\nwith \\(R/L_{n,mid}=%.2f\\)', meta.R_over_Ln_mid);
    chi_solo_formula = [chi_solo_base, rln_mid_text]; % Concatenate parts

    % Define formulas using \( \)
    formula_map('Solomon') = {chi_solo_formula, '\(V_{\mathrm{pinch}}^{(\mathrm{Solomon})} = (-24.2 \pm 3.5)\nu_e^*\)'};
    formula_map('Hahm')    = {chi_solo_formula, '\(V_{\mathrm{pinch}}^{(\mathrm{Hahm})} = -2 \chi_{\varphi} / R_0\)'};
    formula_map('Gurcan')  = {chi_solo_formula, '\(V_{\mathrm{pinch}}^{(\mathrm{G\ddot{u}rcan})} = -(2\chi_{\varphi}/R) \cdot (F+r/R_0)\)'};
    formula_map('Peeters_Rln2') = {chi_solo_formula, '\(V_{\mathrm{pinch}}^{(\mathrm{Peeters})} = (\chi_{\varphi}/R) \cdot (-4 - R/L_n)\), with \(R/L_n=2\)'};
    formula_map('Peeters_Rln_calc') = {chi_solo_formula, '\(V_{\mathrm{pinch}}^{(\mathrm{Peeters})} = (\chi_{\varphi}/R) \cdot (-4 - R/L_n)\)'};
    
    % General formula using $$ remains okay
    general_formula = ['$$\chi_{\varphi}^{\mathrm{eff}} = \chi_{\varphi} \left(1 + \frac{R \cdot V_{\mathrm{pinch}}}{\chi_{\varphi}} ' ...
        '\frac{1}{R/L_{V_{\varphi}}}\right)$$'];
    
    % Define styles with SIMPLE TEXT (no LaTeX) for legend entries
    variant_styles = {
        {'DisplayName', '$R=R_0, R/L_{V_{\varphi}}=%.2f$',    'LineStyle', '--', 'Color', colors(2,:), 'LineWidth', 2.5},
        {'DisplayName', '$R=R_0, R/L_{V_{\varphi}}(r)$',      'LineStyle', '-.', 'Color', colors(3,:), 'LineWidth', 2.0},
        {'DisplayName', '$R=R(r), R/L_{V_{\varphi}}=%.2f$',   'LineStyle', ':',  'Color', colors(4,:), 'LineWidth', 2.5},
        {'DisplayName', '$R=R(r), R/L_{V_{\varphi}}(r)$',     'LineStyle', '--', 'Color', colors(5,:), 'LineWidth', 2.0}
    };
    variant_keys = {'R0_mid', 'R0_profile', 'Rlocal_mid', 'Rlocal_profile'};

    pinch_models = fieldnames(variants);
    pinch_models = pinch_models(~strcmp(pinch_models, 'meta'));
    % --- END OF SECTION 1 ---

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
            'Color', colors(1,:), 'DisplayName', '$\chi_\varphi^{\mathrm{eff}}\ \text{(Exp.)}$');

        % --- 2b. Plot Pure Diffusivity Baseline ---
        plot(ax, r_norm, meta.chi_phi_solo_avg, ':', 'LineWidth', 2, 'Color', 'k', ...
            'DisplayName', '$\chi_\varphi\ \text{(Solomon base)}$');

        % --- Track Max Value for Y-Limit ---
        % Start with max of experimental data (in valid range)
        % We use a range [0.2, 0.95] to avoid edge singularities and core noise
        mask_lims = r_norm >= 0.2 & r_norm <= 0.95;
        
        % Initialize max_val with experiment and base theory
        max_val = max(diffusivity_exp.ci_upper(mask_lims)); 
        max_val = max(max_val, max(meta.chi_phi_solo_avg(mask_lims)));

        % --- 2c. Loop Through and Plot the Four Variants ---
        for j = 1:numel(variant_keys)
            key = variant_keys{j};
            if isfield(model_data, key)
                var_data = model_data.(key);
                style_properties = variant_styles{j};
                
                % Extract the DisplayName template
                display_name_template = '';
                style_args = {};
                is_display_name = false;
                for k = 1:2:numel(style_properties)
                    prop_name = style_properties{k};
                    prop_value = style_properties{k+1};
                    if strcmpi(prop_name, 'DisplayName')
                        display_name_template = prop_value;
                        is_display_name = true;
                    else
                        style_args = [style_args, {prop_name, prop_value}];
                    end
                end
                if ~is_display_name
                    error('DisplayName not found in variant_styles for key %s', key);
                end

                % --- SIMPLE LOGIC FOR LEGEND DISPLAY NAMES ---
                is_mid_variant = any(strcmp(key, {'R0_mid', 'Rlocal_mid'}));
                
                if is_mid_variant
                    if isfield(meta, 'R_over_LVphi_mid') && ~isnan(meta.R_over_LVphi_mid)
                        display_name = sprintf(display_name_template, meta.R_over_LVphi_mid);
                    else
                        display_name = display_name_template;
                    end
                else
                    display_name = display_name_template;
                end

                % Determine Color
                col = [];
                for pp = 1:2:numel(style_properties)
                    if strcmpi(style_properties{pp}, 'Color')
                        col = style_properties{pp+1};
                        break;
                    end
                end
                if isempty(col)
                    col = colors(mod(j-1,size(colors,1))+1,:);
                end

                % Plot
                fill(ax, [r_norm; flipud(r_norm)], [var_data.ci_lower; flipud(var_data.ci_upper)], ...
                     col, 'FaceAlpha', 0.10, 'EdgeColor', 'none', 'HandleVisibility', 'off');
                plot(ax, r_norm, var_data.profile_avg, 'DisplayName', display_name, style_args{:});
                
                % --- Update Max Value for Y-Limit ---
                % Check max of this variant's CI upper bound within mask
                local_max = max(var_data.ci_upper(mask_lims));
                if ~isempty(local_max) && isfinite(local_max)
                    max_val = max(max_val, local_max);
                end
            end
        end

        % --- 2d. Formatting, Legend and Annotations ---
        xlabel('Normalised Radius ($r/a$)', 'Interpreter', 'latex');
        ylabel('Effective Diffusivity, $\chi_{\varphi}^{\mathrm{eff}}$ [m$^2$/s]', 'Interpreter', 'latex');
        title_str = sprintf('Effective Momentum Diffusivity: Experiment vs. %s Scalling Law', strrep(model_name, '_', '\_'));
        title(title_str, 'Interpreter', 'latex');
        
        % Apply Fixed Y-Limit for Consistency
        xlim([0, 1]); 
        ylim([0, 25]); % Fixed limit to accommodate high theoretical values
        
        % Create legend at the top-center
        lgd = legend(ax, 'Location', 'north');
        set(lgd, 'Interpreter', 'latex'); 
        lgd.Box = 'off';
        
        set_publication_style(ax);
        
        % Create the theoretical formula annotation box (keeping layout)
        model_formulas = formula_map(model_name);
        anno_text = { ...
            '\bf{Framework:}', general_formula, ...
            '\bf{Base Diffusivity:}', model_formulas{1}, ...
            '\bf{Pinch Model:}', model_formulas{2} ...
        };
        annotation('textbox', [0.15, 0.44, 0.4, 0.2], 'String', anno_text, ...
            'Interpreter', 'latex', 'FontSize', 16, 'FitBoxToText', 'on', ...
            'BackgroundColor', [1, 1, 1, 0.85], 'EdgeColor', 'k', 'VerticalAlignment', 'top');
        
        hold(ax, 'off');
        
        % --- 2e. Save the Figure ---
        fig_filename = fullfile(save_dir, ['final_comparison_', model_name, '.png']);
        saveas(fig, fig_filename);
        savefig(fig, strrep(fig_filename, '.png', '.fig'));
        fprintf('...plot saved to %s\n', fig_filename);
    end
end