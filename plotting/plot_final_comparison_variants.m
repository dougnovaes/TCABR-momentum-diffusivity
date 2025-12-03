function plot_final_comparison_variants(diffusivity_exp, variants, constants, save_dir)
%PLOT_FINAL_COMPARISON_VARIANTS Generates consolidated comparison plots.
%   Produces two publication-quality figures comparing experimental effective
%   diffusivity with theoretical predictions based on different Pinch mechanisms.
%
%   FIGURES:
%   1. Comparison of Pinch Mechanisms (Solomon, Hahm, Gürcan).
%   2. Sensitivity to Density Gradient (Peeters models).
%
%   VISUALIZATION LOGIC:
%   - Uses fixed variant: Global R0 + Local Velocity Gradient (R0_profile).
%   - Y-axis fixed at [0, 25] m^2/s for consistent scale.
%   - Experimental data and Solomon base: Plotted on FULL range (0 <= r/a <= 1).
%   - Theoretical comparison curves: VISUALLY CROPPED to 0.14 <= r/a <= 0.99
%     to hide non-physical singularities at the core and edge noise, ensuring
%     a clean presentation of the relevant confinement physics.
%
%   Inputs:
%       diffusivity_exp : Struct with experimental profile and CI.
%       variants        : Struct with calculated theoretical variants.
%       constants       : Machine constants.
%       save_dir        : Directory to save output files.

    arguments
        diffusivity_exp (1,1) struct
        variants (1,1) struct
        constants (1,1) struct
        save_dir (1,1) string
    end
    
    fprintf('Generating final consolidated comparison plots...\n');
    
    % --- Setup ---
    r_norm = variants.meta.r_norm;
    colors = get(groot,'DefaultAxesColorOrder');
    meta = variants.meta;
    
    % --- CONFIGURATION: Colors ---
    base_color = 'k'; 
    
    % Define Plotting Limits (Radial Crop) - ONLY for theoretical curves
    R_MIN = 0.14;
    R_MAX = 0.99;
    
    % Helper to compute mask indices for theoretical curves
    get_mask = @(r) (r >= R_MIN) & (r <= R_MAX);
    mask_theo = get_mask(r_norm);
    
    % --- FORMULAS FOR ANNOTATION ---
    % General Framework
    general_formula = '$$\chi_{\phi}^{\mathrm{eff}} = \chi_{\phi} \left(1 + \frac{R \cdot V_{\mathrm{pinch}}}{\chi_{\phi}} \frac{1}{R/L_{V_{\phi}}}\right)$$';
    
    % Base Diffusivity (Solomon)
    chi_solo_formula = sprintf('$\\chi_{\\phi}^{(\\mathrm{Solo})} = (6.09 \\pm 0.72)\\nu_e^* + (0.157 \\pm 0.072)R/L_n$,\nwith $R/L_{n,mid}=%.2f$', meta.R_over_Ln_mid);

    % --- Common Plotting Helper ---
    function plot_base_layers(ax)
        % 1. Experimental Data (Blue band and line) - FULL RANGE
        if isfield(diffusivity_exp, 'r_fine')
            r_exp_norm = diffusivity_exp.r_fine / constants.machine.a;
        else
            r_exp_norm = r_norm; 
        end
        % Experimental CI
        fill(ax, [r_exp_norm; flipud(r_exp_norm)], [diffusivity_exp.ci_lower; flipud(diffusivity_exp.ci_upper)], ...
            colors(1,:), 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'DisplayName', '95\% CI (Exp.)');
        % Experimental Mean
        plot(ax, r_exp_norm, diffusivity_exp.profile_avg, '-', 'LineWidth', 3.5, ...
            'Color', colors(1,:), 'DisplayName', '$\chi_\phi^{\mathrm{eff}}\ \text{(Exp.)}$');
        
        % 2. Base Diffusivity (Solomon) - FULL RANGE with Confidence Band
        if isfield(meta, 'chi_phi_solo_lower')
            fill(ax, [r_norm; flipud(r_norm)], [meta.chi_phi_solo_lower; flipud(meta.chi_phi_solo_upper)], ...
                 [0.6 0.6 0.6], 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        end
        % Base Mean Line
        plot(ax, r_norm, meta.chi_phi_solo_avg, ':', 'LineWidth', 2, 'Color', base_color, ...
            'DisplayName', '$\chi_\phi\ \text{(Solomon base)}$');
    end

    % =====================================================================
    % FIGURE 1: COMPARISON OF PINCH MECHANISMS (Solomon, Hahm, Gürcan)
    % =====================================================================
    fig1 = figure('Name', 'Comparison: Pinch Mechanisms', 'WindowStyle', 'docked');
    ax1 = gca; hold(ax1, 'on');
    
    plot_base_layers(ax1);
    
    % Define Models to Plot
    mech_models = {'Solomon', 'Hahm', 'Gurcan'};
    mech_colors = [0.9290 0.6940 0.1250;  % Yellow (Solomon)
                   0.4940 0.1840 0.5560;  % Purple (Hahm)
                   0.4660 0.6740 0.1880]; % Green (Gürcan)
    mech_linewidths = [2.5, 4.0, 2.5];
    
    % Loop and Plot Variants (CROPPED)
    for i = 1:numel(mech_models)
        name = mech_models{i};
        if isfield(variants, name)
            % Use "R0_profile" variant (Global R, Local Grads)
            data = variants.(name).R0_profile; 
            col = mech_colors(i,:);
            
            % Apply Mask - CROPPED theoretical curves
            r_t = r_norm(mask_theo);
            prof_t = data.profile_avg(mask_theo);
            ci_lo_t = data.ci_lower(mask_theo);
            ci_up_t = data.ci_upper(mask_theo);
            
            % Plot CI (faint)
            fill(ax1, [r_t; flipud(r_t)], [ci_lo_t; flipud(ci_up_t)], ...
                 col, 'FaceAlpha', 0.08, 'EdgeColor', 'none', 'HandleVisibility', 'off');
            % Plot Mean Line
            plot(ax1, r_t, prof_t, '--', 'LineWidth', mech_linewidths(i), 'Color', col, ...
                 'DisplayName', sprintf('%s Model', name));
        end
    end
    
    % Formatting
    xlabel(ax1, 'Normalised Radius ($r/a$)', 'Interpreter', 'latex');
    ylabel(ax1, 'Effective Diffusivity [m$^2$/s]', 'Interpreter', 'latex');
    title(ax1, 'Comparison of Pinch Mechanisms (TEP vs Coriolis)', 'Interpreter', 'latex');
    
    xlim(ax1, [0, 1]); ylim(ax1, [0, 25]);
    
    legend(ax1, 'Location', 'southeast', 'Interpreter', 'latex', 'Box', 'off', 'NumColumns', 1);
    leg1 = legend(ax1, 'Interpreter', 'latex', 'Box', 'off');
    leg1.Position = [0.45, 0.6, 0.30, 0.35]; % [left, bottom, width, height]
    set_publication_style(ax1);
    
    % Annotation Box (UPDATED FORMULAS)
    anno_text1 = { ...
        '\bf{Framework:}', general_formula, ...
        '\bf{Base Diffusivity:}', chi_solo_formula, ...
        '\bf{Pinch Models:}', ...
        '$V_{\mathrm{pinch}}^{(\mathrm{Solo})} = (-24.2 \pm 3.5)\nu_e^*$', ...
        '$V_{\mathrm{pinch}}^{(\mathrm{Hahm})} = -2 (\chi_{\phi} / R_0) \cdot F$', ...
        '$V_{\mathrm{pinch}}^{(\mathrm{G\ddot{u}rc})} = \frac{\chi_{\phi}}{R} [-2(F + r/R_0) + R/L_n]$'
    };
    annotation('textbox', [0.14, 0.55, 0.38, 0.35], 'String', anno_text1, ...
        'Interpreter', 'latex', 'FontSize', 15, 'FitBoxToText', 'on', ...
        'BackgroundColor', [1, 1, 1, 0.85], 'EdgeColor', 'k', 'VerticalAlignment', 'top');
    
    saveas(fig1, fullfile(save_dir, 'final_comparison_mechanisms.png'));
    savefig(fig1, fullfile(save_dir, 'final_comparison_mechanisms.fig'));


    % =====================================================================
    % FIGURE 2: SENSITIVITY TO DENSITY GRADIENT (Peeters)
    % =====================================================================
    fig2 = figure('Name', 'Comparison: Peeters Sensitivity', 'WindowStyle', 'docked');
    ax2 = gca; hold(ax2, 'on');
    
    plot_base_layers(ax2);
    
    % Prepare common masking variables for theoretical curves (CROPPED)
    r_t = r_norm(mask_theo);
    
    % Peeters R/Ln = 2 (Standard) - CROPPED
    if isfield(variants, 'Peeters_Rln2')
        data = variants.Peeters_Rln2.R0_profile;
        col = [0.8500 0.3250 0.0980]; % Red-Orange
        
        prof = data.profile_avg(mask_theo);
        ci_lo = data.ci_lower(mask_theo);
        ci_up = data.ci_upper(mask_theo);
        
        fill(ax2, [r_t; flipud(r_t)], [ci_lo; flipud(ci_up)], ...
                 col, 'FaceAlpha', 0.08, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        plot(ax2, r_t, prof, '-.', 'LineWidth', 2.5, 'Color', col, ...
             'DisplayName', 'Peeters ($R/L_n = 2$)');
    end
    
    % Peeters R/Ln Calculated (Experimental ~ 6) - CROPPED
    if isfield(variants, 'Peeters_Rln_calc')
        data = variants.Peeters_Rln_calc.R0_profile;
        col = [0.6350 0.0780 0.1840]; % Dark Red
        
        prof = data.profile_avg(mask_theo);
        ci_lo = data.ci_lower(mask_theo);
        ci_up = data.ci_upper(mask_theo);
        
        fill(ax2, [r_t; flipud(r_t)], [ci_lo; flipud(ci_up)], ...
                 col, 'FaceAlpha', 0.08, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        plot(ax2, r_t, prof, '--', 'LineWidth', 2.5, 'Color', col, ...
             'DisplayName', sprintf('Peeters ($R/L_n \\approx %.1f$)', meta.R_over_Ln_mid));
    end
    
    % Formatting
    xlabel(ax2, 'Normalised Radius ($r/a$)', 'Interpreter', 'latex');
    ylabel(ax2, 'Effective Diffusivity [m$^2$/s]', 'Interpreter', 'latex');
    title(ax2, 'Peeters Model: Sensitivity to Density Gradient', 'Interpreter', 'latex');
    
    xlim(ax2, [0, 1]); ylim(ax2, [0, 25]);
    
    % legend(ax2, 'Location', 'northeast', 'Interpreter', 'latex', 'Box', 'off');
    leg2 = legend(ax2, 'Interpreter', 'latex', 'Box', 'off');
    leg2.Position = [0.45, 0.7, 0.25, 0.15]; % [left, bottom, width, height]
    set_publication_style(ax2);
    
    % Annotation Box (PEETERS)
    anno_text2 = { ...
        '\bf{Framework:}', general_formula, ...
        '\bf{Base Diffusivity:}', chi_solo_formula, ...
        '\bf{Pinch Model (Peeters):}', ...
        '$V_{\mathrm{pinch}} = (\chi_{\phi}/R) \cdot (-4 - R/L_n)$', ...
        'Testing sensitivity to $R/L_n$ choice.'
    };
    annotation('textbox', [0.14, 0.65, 0.38, 0.25], 'String', anno_text2, ...
        'Interpreter', 'latex', 'FontSize', 15, 'FitBoxToText', 'on', ...
        'BackgroundColor', [1, 1, 1, 0.85], 'EdgeColor', 'k', 'VerticalAlignment', 'top');

    saveas(fig2, fullfile(save_dir, 'final_comparison_peeters.png'));
    savefig(fig2, fullfile(save_dir, 'final_comparison_peeters.fig'));
    
    fprintf('...consolidated plots saved.\n');
end