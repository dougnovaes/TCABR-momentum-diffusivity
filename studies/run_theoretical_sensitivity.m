function run_theoretical_sensitivity(velocity_results, derived_profiles, constants, r_fine, results_dir)
%RUN_THEORETICAL_SENSITIVITY Explores parameter sensitivity in theoretical models.
%   Focuses on:
%   1. Gürcan Model: Sensitivity to poloidal asymmetry factor 'F'.
%   2. Peeters Model: Sensitivity to density gradient scale length R/Ln.

    arguments
        velocity_results (1,1) struct
        derived_profiles (1,1) struct
        constants (1,1) struct
        r_fine (:,1) double
        results_dir (1,1) string
    end

    fprintf('\n=== Starting Theoretical Sensitivity Study ===\n');
    
    % --- Setup ---
    a = constants.machine.a;
    R0 = constants.machine.R0;
    r_norm = r_fine / a;
    R_coord = R0 * (1 + r_fine / R0);
    
    % Base Inputs
    % Use Solomon collisionality if available (Act 3 physics)
    if isfield(derived_profiles, 'nu_star_e_solomon')
        nu_star_e = derived_profiles.nu_star_e_solomon(:);
    else
        nu_star_e = derived_profiles.nu_star_e(:);
    end
    
    % Gradients
    v_phi = velocity_results.poly_fit_avg;
    grad_v = gradient(v_phi, r_fine);
    v_safe = v_phi; v_safe(abs(v_safe)<1e-9) = 1e-9;
    R_over_LVphi = -R0 ./ v_safe .* grad_v;
    
    R_over_Ln_local = derived_profiles.R_over_Ln;
    [~, idx_mid] = min(abs(r_fine - a*0.5));
    R_over_Ln_mid = derived_profiles.R_over_Ln(idx_mid);

    % Base Diffusivity (Solomon) - Fixed for comparisons
    A = 6.09; B = 0.157;
    chi_phi_base = A * nu_star_e + B * R_over_Ln_mid; % Using mid for baseline stability
    
    % Output dir
    sens_dir = fullfile(results_dir, 'sensitivity_study');
    if ~exist(sens_dir, 'dir'), mkdir(sens_dir); end

    % =====================================================================
    % STUDY 1: GÜRCAN MODEL - SENSITIVITY TO 'F'
    % V_pinch = -(2*chi/R) * (F + r/R0)
    % =====================================================================
    fprintf('>> Analyzing Gürcan Model sensitivity to F...\n');
    F_values = [0.0, 0.5, 1.0];
    colors_F = lines(length(F_values));
    
    f1 = figure('Name', 'Sensitivity: Gurcan F', 'WindowStyle', 'docked');
    ax1 = gca; hold on; grid on; box on;
    
    % Plot Base Chi (Limit of no pinch effect)
    plot(ax1, r_norm, chi_phi_base, 'k--', 'LineWidth', 1.5, 'DisplayName', '\chi_{\phi} (Pure Diffusive)');
    
    for i = 1:length(F_values)
        F = F_values(i);
        % Calculate Pinch
        v_pinch = -(2 * chi_phi_base ./ R_coord) .* (F + r_fine./R0);
        
        % Calculate Effective Diffusivity
        % chi_eff = chi * (1 + R*Vp/chi * 1/(R/Lv))
        term_pinch = (R0 * v_pinch ./ chi_phi_base) ./ R_over_LVphi;
        chi_eff = chi_phi_base .* (1 + term_pinch);
        
        % Handle singularities for plotting
        chi_eff(abs(chi_eff) > 50) = NaN; 
        
        plot(ax1, r_norm, chi_eff, 'LineWidth', 2, 'Color', colors_F(i,:), ...
            'DisplayName', sprintf('F = %.1f', F));
    end
    
    title(ax1, 'Gürcan Model Sensitivity to Asymmetry Factor F');
    xlabel(ax1, 'r/a'); ylabel(ax1, '\chi_{\phi, eff} [m^2/s]');
    legend(ax1, 'Location', 'best');
    ylim(ax1, [0, 30]); xlim(ax1, [0.2, 0.9]);
    saveas(f1, fullfile(sens_dir, 'sensitivity_gurcan_F.png'));


    % =====================================================================
    % STUDY 2: PEETERS MODEL - SENSITIVITY TO R/Ln
    % V_pinch = (chi/R) * (-4 - R/Ln)
    % =====================================================================
    fprintf('>> Analyzing Peeters Model sensitivity to R/Ln...\n');
    
    f2 = figure('Name', 'Sensitivity: Peeters R/Ln', 'WindowStyle', 'docked');
    ax2 = gca; hold on; grid on; box on;
    
    plot(ax2, r_norm, chi_phi_base, 'k--', 'LineWidth', 1.5, 'DisplayName', '\chi_{\phi} (Pure Diffusive)');
    
    % Case A: Constant R/Ln = 2 (Standard)
    vp_2 = (chi_phi_base ./ R_coord) .* (-4 - 2.0);
    chi_eff_2 = chi_phi_base .* (1 + (R0 * vp_2 ./ chi_phi_base) ./ R_over_LVphi);
    chi_eff_2(abs(chi_eff_2)>50) = NaN;
    plot(ax2, r_norm, chi_eff_2, 'b-', 'LineWidth', 2, 'DisplayName', 'R/L_n = 2.0 (Fixed)');
    
    % Case B: Constant R/Ln = 6 (TCABR Mid)
    vp_6 = (chi_phi_base ./ R_coord) .* (-4 - 6.0);
    chi_eff_6 = chi_phi_base .* (1 + (R0 * vp_6 ./ chi_phi_base) ./ R_over_LVphi);
    chi_eff_6(abs(chi_eff_6)>50) = NaN;
    plot(ax2, r_norm, chi_eff_6, 'r-', 'LineWidth', 2, 'DisplayName', 'R/L_n = 6.0 (Fixed)');
    
    % Case C: Local Profile R/Ln(r)
    vp_prof = (chi_phi_base ./ R_coord) .* (-4 - R_over_Ln_local);
    chi_eff_prof = chi_phi_base .* (1 + (R0 * vp_prof ./ chi_phi_base) ./ R_over_LVphi);
    chi_eff_prof(abs(chi_eff_prof)>50) = NaN;
    plot(ax2, r_norm, chi_eff_prof, 'g-', 'LineWidth', 2, 'DisplayName', 'R/L_n(r) (Local Profile)');
    
    title(ax2, 'Peeters Model Sensitivity to Density Gradient R/L_n');
    xlabel(ax2, 'r/a'); ylabel(ax2, '\chi_{\phi, eff} [m^2/s]');
    legend(ax2, 'Location', 'best');
    ylim(ax2, [0, 30]); xlim(ax2, [0.2, 0.9]);
    saveas(f2, fullfile(sens_dir, 'sensitivity_peeters_RLn.png'));
    
    fprintf('Sensitivity plots saved to: %s\n', sens_dir);
end