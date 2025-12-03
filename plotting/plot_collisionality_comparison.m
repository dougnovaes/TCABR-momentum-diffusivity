function plot_collisionality_comparison(derived_profiles, constants, velocity_results, diffusivity_exp, save_dir)
%PLOT_COLLISIONALITY_COMPARISON Visualizes the impact of collisionality definitions.
%
%   PURPOSE:
%       Generates a 3-panel plot to demonstrate the sensitivity of the 
%       theoretical diffusivity to the choice of collisionality definition
%       (Wesson/Bounce vs. Solomon/Drift).
%
%   PANELS:
%       (a) Collisionality Profiles (\nu*).
%       (b) Base Diffusivity Profiles (\chi_phi).
%       (c) Effective Diffusivity Profiles (\chi_eff) compared to experiment.
%
%   INPUTS:
%       derived_profiles : Struct containing collisionality profiles.
%       constants        : Project constants.
%       velocity_results : Struct with Vphi fit (for R/Lv).
%       diffusivity_exp  : Experimental diffusivity results.
%       save_dir         : Output directory.

    arguments
        derived_profiles (1,1) struct
        constants (1,1) struct
        velocity_results (1,1) struct
        diffusivity_exp (1,1) struct
        save_dir (1,1) string
    end

    fprintf('Generating collisionality definition comparison plot...\n');

    % --- 1. Check Availability ---
    if ~isfield(derived_profiles, 'nu_star_e_solomon')
        warning('Solomon collisionality not found. Skipping plot.');
        return;
    end

    % --- 2. Recalculate Simple Variants for Display ---
    a = constants.machine.a;
    R0 = constants.machine.R0;
    [~, idx_mid] = min(abs(derived_profiles.r_fine - a * 0.5));
    
    if isfield(derived_profiles, 'R_over_Ln')
        R_over_Ln_mid = derived_profiles.R_over_Ln(idx_mid);
    else
        R_over_Ln_mid = derived_profiles.gradients.R_over_Ln(idx_mid);
    end
    
    v_phi_smooth = smoothdata(velocity_results.poly_fit_avg, 'movmean', 15);
    grad_Vphi = gradient(v_phi_smooth, derived_profiles.r_fine);
    v_phi_safe = v_phi_smooth; v_phi_safe(abs(v_phi_safe) < 1e-3) = 1e-3;
    R_over_LVphi_mid = (-R0 ./ v_phi_safe(idx_mid)) .* grad_Vphi(idx_mid);

    A = 6.09; B = 0.157; C = -24.2;
    calc_chi_eff = @(nu) (A * nu + B * R_over_Ln_mid) .* ...
                         (1 + (R0 * C * nu) ./ (A * nu + B * R_over_Ln_mid) / R_over_LVphi_mid);

    % Data Extraction
    if isfield(derived_profiles, 'nu_star_e_wesson')
        nu_wesson = derived_profiles.nu_star_e_wesson;
    else
        nu_wesson = derived_profiles.nu_star_e;
    end
    
    % For this specific plot, we want to show the LOCAL profile of Solomon 
    % to demonstrate why it explodes, contrasting with the global approach used later.
    if isfield(derived_profiles, 'nu_star_e_solomon_local')
        nu_solomon = derived_profiles.nu_star_e_solomon_local;
    else
        nu_solomon = derived_profiles.nu_star_e_solomon;
    end
    
    chi_phi_wesson = A * nu_wesson + B * R_over_Ln_mid;
    chi_phi_solomon = A * nu_solomon + B * R_over_Ln_mid;
    
    chi_eff_wesson = calc_chi_eff(nu_wesson);
    chi_eff_solomon = calc_chi_eff(nu_solomon);

    % --- 3. Plotting ---
    fig = figure('Name', 'Collisionality Sensitivity Analysis', 'WindowStyle', 'docked');
    r_norm = derived_profiles.r_fine / a;
    
    % Panel (a)
    ax1 = subplot(1, 3, 1); hold(ax1, 'on'); grid(ax1, 'on'); box(ax1, 'on');
    plot(ax1, r_norm, nu_wesson, 'LineWidth', 2.5, 'DisplayName', 'Wesson (Bounce)');
    plot(ax1, r_norm, nu_solomon, 'LineWidth', 2.5, 'DisplayName', 'Solomon (Drift)');
    xlabel(ax1, 'r/a'); ylabel(ax1, '$\nu_{e}^*$'); 
    title(ax1, '(a) Electron Collisionality');
    legend(ax1, 'Location', 'best'); xlim(ax1, [0, 1]); set_publication_style(ax1);

    % Panel (b)
    ax2 = subplot(1, 3, 2); hold(ax2, 'on'); grid(ax2, 'on'); box(ax2, 'on');
    plot(ax2, r_norm, chi_phi_wesson, 'LineWidth', 2.5, 'DisplayName', 'Using \nu_{Wesson}');
    plot(ax2, r_norm, chi_phi_solomon, 'LineWidth', 2.5, 'DisplayName', 'Using \nu_{Solomon}');
    xlabel(ax2, 'r/a'); ylabel(ax2, '$\chi_{\phi} [m^2/s]$'); 
    title(ax2, '(b) Base Diffusivity $\chi_{\phi}$');
    legend(ax2, 'Location', 'best'); xlim(ax2, [0, 1]); set_publication_style(ax2);
    
    % Panel (c)
    ax3 = subplot(1, 3, 3); hold(ax3, 'on'); grid(ax3, 'on'); box(ax3, 'on');
    fill(ax3, [r_norm; flipud(r_norm)], [diffusivity_exp.ci_lower; flipud(diffusivity_exp.ci_upper)], ...
        [0.8 0.8 0.8], 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'DisplayName', 'Exp. 95\% CI');
    plot(ax3, r_norm, diffusivity_exp.profile_avg, 'k-', 'LineWidth', 2, 'DisplayName', 'Exp. Mean');
    plot(ax3, r_norm, chi_eff_wesson, '--', 'LineWidth', 2, 'DisplayName', 'Theory (\nu_{Wesson})');
    plot(ax3, r_norm, chi_eff_solomon, '-.', 'LineWidth', 2, 'DisplayName', 'Theory (\nu_{Solomon})');
    xlabel(ax3, 'r/a'); ylabel(ax3, '$\chi_{\phi, \mathrm{eff}} [m^2/s]$'); 
    title(ax3, '(c) Effective Diffusivity Sensitivity');
    legend(ax3, 'Location', 'best'); xlim(ax3, [0, 1]); ylim(ax3, [0, 60]); 
    set_publication_style(ax3);

    if ~exist(save_dir, 'dir'), mkdir(save_dir); end
    saveas(fig, fullfile(save_dir, 'analysis_collisionality_comparison.png'));
    savefig(fig, fullfile(save_dir, 'analysis_collisionality_comparison.fig'));
end