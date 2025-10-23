function plot_collisionality_comparison(derived_profiles, constants, velocity_results, diffusivity_exp, save_dir)
%PLOT_COLLISIONALITY_COMPARISON Visualizes the impact of different collisionality definitions.
%   This function generates a 3-panel plot to explicitly compare and quantify
%   the effect of using the Wesson-based vs. the Solomon/Maslov-based
%   electron collisionality definition (`nu_star_e`).
%
%   Panels:
%   (a) Direct comparison of the two `nu_star_e` profiles.
%   (b) The resulting impact on the pure momentum diffusivity, `chi_phi`.
%   (c) The final impact on the effective momentum diffusivity, `chi_phi_eff`,
%       shown against the experimental result.

    arguments
        derived_profiles (1,1) struct
        constants (1,1) struct
        velocity_results (1,1) struct
        diffusivity_exp (1,1) struct
        save_dir (1,1) string
    end

    fprintf('Generating collisionality definition comparison plot...\n');

    % --- 1. Recalculate Theoretical Variants for Both Collisionalities ---
    % We will calculate a simple theoretical case (Solomon pinch, R0, mid-grads)
    % for each definition of nu_star_e to show the impact.
    
    % Extract common parameters
    a = constants.machine.a;
    R0 = constants.machine.R0;
    [~, idx_mid] = min(abs(derived_profiles.r_fine - a * 0.5));
    R_over_Ln_mid = derived_profiles.gradients.R_over_Ln(idx_mid);
    
    % Use the robustly calculated R/L_Vphi from the variants metadata
    % (Requires Stage 7b to have been run at least once)
    if isfield(derived_profiles, 'variants') && isfield(derived_profiles.variants.meta, 'R_over_LVphi_mid')
        R_over_LVphi_mid = derived_profiles.variants.meta.R_over_LVphi_mid;
    else % Fallback calculation if variants are not available
        v_phi_smooth = smoothdata(velocity_results.poly_fit_avg, 'movmean', 15);
        grad_Vphi = gradient(v_phi_smooth, derived_profiles.r_fine);
        v_phi_safe = v_phi_smooth; v_phi_safe(abs(v_phi_safe) < 1e-3) = 1e-3;
        R_over_LVphi_profile = -R0 ./ v_phi_safe .* grad_Vphi;
        R_over_LVphi_mid = R_over_LVphi_profile(idx_mid);
    end

    % Solomon coefficients
    A = 6.09; B = 0.157; C = -24.2;

    % Function to calculate chi_eff from a given nu_star_e
    calculate_chi = @(nu_star_e) ...
        (A * nu_star_e + B * R_over_Ln_mid) .* ...
        (1 + (R0 * C * nu_star_e) ./ (A * nu_star_e + B * R_over_Ln_mid) / R_over_LVphi_mid);

    % Calculate for both cases
    chi_phi_wesson = A * derived_profiles.collisionality.nu_star_e_wesson + B * R_over_Ln_mid;
    chi_phi_solomon = A * derived_profiles.collisionality.nu_star_e_solomon + B * R_over_Ln_mid;
    
    chi_eff_wesson = calculate_chi(derived_profiles.collisionality.nu_star_e_wesson);
    chi_eff_solomon = calculate_chi(derived_profiles.collisionality.nu_star_e_solomon);

    % --- 2. Create Figure and Subplots ---
    fig = figure('Name', 'Collisionality Definition Comparison', 'Position', [100, 100, 1800, 600]);
    r_norm = derived_profiles.r_fine / a;
    
    % Panel (a): Collisionality Profiles
    ax1 = subplot(1, 3, 1); hold(ax1, 'on'); grid(ax1, 'on'); box(ax1, 'on');
    plot(ax1, r_norm, derived_profiles.collisionality.nu_star_e_wesson, 'LineWidth', 3, 'DisplayName', 'Wesson (bounce freq.)');
    plot(ax1, r_norm, derived_profiles.collisionality.nu_star_e_solomon, 'LineWidth', 3, 'DisplayName', 'Solomon/Maslov (drift freq.)');
    xlabel(ax1, 'Normalised Radius ($r/a$)'); ylabel(ax1, 'Electron Collisionality, $\nu_{e}^*$');
    title(ax1, '(a) Collisionality Definitions');
    legend(ax1, 'show', 'Location', 'northwest');
    set_publication_style(ax1);

    % Panel (b): Impact on Pure Diffusivity
    ax2 = subplot(1, 3, 2); hold(ax2, 'on'); grid(ax2, 'on'); box(ax2, 'on');
    plot(ax2, r_norm, chi_phi_wesson, 'LineWidth', 3, 'DisplayName', '$\chi_\varphi$ using Wesson $\nu_e^*$');
    plot(ax2, r_norm, chi_phi_solomon, 'LineWidth', 3, 'DisplayName', '$\chi_\varphi$ using Solomon $\nu_e^*$');
    xlabel(ax2, 'Normalised Radius ($r/a$)'); ylabel(ax2, 'Pure Diffusivity, $\chi_\varphi$ [m$^2$/s]');
    title(ax2, '(b) Impact on Pure Diffusivity');
    legend(ax2, 'show', 'Location', 'northwest');
    set_publication_style(ax2);
    
    % Panel (c): Impact on Effective Diffusivity
    ax3 = subplot(1, 3, 3); hold(ax3, 'on'); grid(ax3, 'on'); box(ax3, 'on');
    plot(ax3, r_norm, diffusivity_exp.profile_avg, 'k-', 'LineWidth', 3.5, 'DisplayName', 'Experimental');
    plot(ax3, r_norm, chi_eff_wesson, 'LineWidth', 3, 'DisplayName', 'Theory using Wesson $\nu_e^*$');
    plot(ax3, r_norm, chi_eff_solomon, 'LineWidth', 3, 'DisplayName', 'Theory using Solomon $\nu_e^*$');
    xlabel(ax3, 'Normalised Radius ($r/a$)'); ylabel(ax3, 'Effective Diffusivity, $\chi_{\phi, \mathrm{eff}}$ [m$^2$/s]');
    title(ax3, '(c) Impact on Effective Diffusivity');
    legend(ax3, 'show', 'Location', 'northwest');
    set_publication_style(ax3);

    % --- 3. Save Figure ---
    if ~exist(save_dir, 'dir'), mkdir(save_dir); end
    fig_filename = fullfile(save_dir, 'analysis_collisionality_comparison');
    saveas(fig, [fig_filename, '.png']);
    savefig(fig, [fig_filename, '.fig']);
    fprintf('...collisionality comparison plot saved to %s.(png/fig)\n', fig_filename);
end