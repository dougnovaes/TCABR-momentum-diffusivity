function plot_collisionality_comparison(derived_profiles, constants, velocity_results, diffusivity_exp, save_dir)
%PLOT_COLLISIONALITY_COMPARISON Visualizes the impact of different collisionality definitions.
%   Generates a 3-panel plot comparing Wesson (bounce freq) vs Solomon (curvature drift)
%   collisionality and their propagation into chi_phi and chi_eff.

    arguments
        derived_profiles (1,1) struct
        constants (1,1) struct
        velocity_results (1,1) struct
        diffusivity_exp (1,1) struct
        save_dir (1,1) string
    end

    fprintf('Generating collisionality definition comparison plot...\n');

    % --- 1. Check availability of Solomon collisionality ---
    if ~isfield(derived_profiles, 'nu_star_e_solomon')
        warning('Solomon collisionality (nu_star_e_solomon) not found. Skipping plot.');
        return;
    end

    % --- 2. Recalculate Simple Variants ---
    a = constants.machine.a;
    R0 = constants.machine.R0;
    [~, idx_mid] = min(abs(derived_profiles.r_fine - a * 0.5));
    
    % Get gradients (Robust access)
    if isfield(derived_profiles, 'R_over_Ln')
        R_over_Ln_mid = derived_profiles.R_over_Ln(idx_mid);
    else
        R_over_Ln_mid = derived_profiles.gradients.R_over_Ln(idx_mid);
    end
    
    % Re-estimate R/L_Vphi
    v_phi_smooth = smoothdata(velocity_results.poly_fit_avg, 'movmean', 15);
    grad_Vphi = gradient(v_phi_smooth, derived_profiles.r_fine);
    v_phi_safe = v_phi_smooth; v_phi_safe(abs(v_phi_safe) < 1e-3) = 1e-3;
    R_over_LVphi_profile = -R0 ./ v_phi_safe .* grad_Vphi;
    R_over_LVphi_mid = R_over_LVphi_profile(idx_mid);

    % Solomon coefficients
    A = 6.09; B = 0.157; C = -24.2;

    % Helper to calc chi_eff
    calc_chi_eff = @(nu) (A * nu + B * R_over_Ln_mid) .* ...
                         (1 + (R0 * C * nu) ./ (A * nu + B * R_over_Ln_mid) / R_over_LVphi_mid);

    % Data
    if isfield(derived_profiles, 'nu_star_e_wesson')
        nu_wesson = derived_profiles.nu_star_e_wesson;
    else
        nu_wesson = derived_profiles.nu_star_e;
    end
    nu_solomon = derived_profiles.nu_star_e_solomon;
    
    chi_phi_wesson = A * nu_wesson + B * R_over_Ln_mid;
    chi_phi_solomon = A * nu_solomon + B * R_over_Ln_mid;
    
    chi_eff_wesson = calc_chi_eff(nu_wesson);
    chi_eff_solomon = calc_chi_eff(nu_solomon);

    % --- 3. Plotting ---
    fig = figure('Name', 'Collisionality Sensitivity Analysis', 'WindowStyle', 'docked');
    r_norm = derived_profiles.r_fine / a;
    
    % Panel (a): Collisionality Profiles
    ax1 = subplot(1, 3, 1); hold(ax1, 'on'); grid(ax1, 'on'); box(ax1, 'on');
    plot(ax1, r_norm, nu_wesson, 'LineWidth', 2.5, 'DisplayName', 'Wesson (Bounce)');
    plot(ax1, r_norm, nu_solomon, 'LineWidth', 2.5, 'DisplayName', 'Solomon (Drift)');
    
    % FIX: Use single backslash for LaTeX commands in standard strings
    xlabel(ax1, 'r/a'); ylabel(ax1, '$\nu_{e}^*$'); 
    title(ax1, '(a) Electron Collisionality');
    legend(ax1, 'Location', 'best');
    xlim(ax1, [0, 1]); set_publication_style(ax1);

    % Panel (b): Pure Diffusivity
    ax2 = subplot(1, 3, 2); hold(ax2, 'on'); grid(ax2, 'on'); box(ax2, 'on');
    plot(ax2, r_norm, chi_phi_wesson, 'LineWidth', 2.5, 'DisplayName', 'Using \nu_{Wesson}');
    plot(ax2, r_norm, chi_phi_solomon, 'LineWidth', 2.5, 'DisplayName', 'Using \nu_{Solomon}');
    
    % FIX: Single backslash
    xlabel(ax2, 'r/a'); ylabel(ax2, '$\chi_{\varphi} [\mathrm{m^2/s}]$');
    title(ax2, '(b) Base Diffusivity $\chi_{\varphi}$');
    legend(ax2, 'Location', 'best');
    xlim(ax2, [0, 1]); set_publication_style(ax2);
    % 'Effective Diffusivity, $\chi_{\varphi}^{\mathrm{eff}}$ [m$^2$/s]', 'Interpreter', 'latex'
    % Panel (c): Effective Diffusivity vs Exp
    ax3 = subplot(1, 3, 3); hold(ax3, 'on'); grid(ax3, 'on'); box(ax3, 'on');
    
    % Exp data with CI
    fill(ax3, [r_norm; flipud(r_norm)], [diffusivity_exp.ci_lower; flipud(diffusivity_exp.ci_upper)], ...
        [0.8 0.8 0.8], 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'DisplayName', 'Exp. 95% CI');
    plot(ax3, r_norm, diffusivity_exp.profile_avg, 'k-', 'LineWidth', 2, 'DisplayName', 'Exp. Mean');
    
    plot(ax3, r_norm, chi_eff_wesson, '--', 'LineWidth', 2, 'DisplayName', 'Theory (\nu_{Wesson})');
    plot(ax3, r_norm, chi_eff_solomon, '-.', 'LineWidth', 2, 'DisplayName', 'Theory (\nu_{Solomon})');
    
    % FIX: Single backslash
    xlabel(ax3, 'r/a'); ylabel(ax3, '$\chi_{\varphi, \mathrm{eff}} [\mathrm{m^2/s}]$');
    title(ax3, '(c) Effective Diffusivity Sensitivity');
    legend(ax3, 'Location', 'best');
    
    % Adjust ylim to accommodate the higher Solomon values
    xlim(ax3, [0, 1]); ylim(ax3, [0, 25]); 
    set_publication_style(ax3);

    % Save
    if ~exist(save_dir, 'dir'), mkdir(save_dir); end
    saveas(fig, fullfile(save_dir, 'analysis_collisionality_comparison.png'));
    savefig(fig, fullfile(save_dir, 'analysis_collisionality_comparison.fig'));
    fprintf('...collisionality comparison plot saved.\n');
end