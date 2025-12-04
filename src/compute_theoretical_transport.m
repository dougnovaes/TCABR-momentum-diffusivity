function [theory_results] = compute_theoretical_transport(velocity_results, derived_profiles, constants, r_fine)
%COMPUTE_THEORETICAL_TRANSPORT Computes theoretical momentum transport coefficients.
%
%   PURPOSE:
%       Calculates the "Pure" Turbulent Diffusivity (Solomon) and Effective
%       Diffusivity (chi_eff) profiles for various Pinch models (Hahm, Gürcan, Peeters).
%
%   METHODOLOGY (HYBRID GLOBAL/LOCAL):
%       1. MAGNITUDE (Global): Uses GLOBAL representative parameters (Volume-averaged
%          <Te>, <ne>, resulting in Maslov's nu*) to determine the scalar baseline
%          levels for diffusivity and pinch velocity. This prevents edge singularities.
%       2. SHAPE (Local): Uses the LOCAL experimental velocity gradient (R/L_Vphi)
%          to modulate the effective diffusivity, capturing the transport topology.
%       3. UNCERTAINTY: Performs a Monte Carlo propagation of uncertainties in
%          the empirical coefficients (Solomon et al., 2010).
%
%   INPUTS:
%       velocity_results : Contains V_phi polynomial fit (for local gradients).
%       derived_profiles : Contains global parameters (nu*_global, R/Ln_maslov).
%
%   OUTPUTS:
%       theory_results   : Struct with mean profiles, CI bands, and scalar diagnostics.

    arguments
        velocity_results (1,1) struct
        derived_profiles (1,1) struct
        constants (1,1) struct
        r_fine (:,1) {mustBeNumeric, mustBeReal, mustBeFinite}
    end

    fprintf('Calculating theoretical transport (Hybrid Global/Local approach)...\n');
    t_start = tic;

    % =========================================================================
    % 1. GEOMETRY & LOCAL GRADIENTS (SHAPE)
    % =========================================================================
    a  = constants.machine.a;
    R0 = constants.machine.R0;
    r_norm = r_fine / a;
    R_coord = R0 * (1 + r_fine / R0); % Local Major Radius
    epsilon_global = constants.machine.epsilon_aspect_ratio;

    % Velocity Gradient Profile (R/L_Vphi)
    v_phi = velocity_results.poly_fit_avg; 
    grad_Vphi = gradient(v_phi, r_fine);
    v_safe = v_phi; 
    v_safe(abs(v_safe) < 1e-9) = 1e-9; % Avoid division by zero
    
    % This profile defines the "shape" of the effective diffusivity
    R_over_LVphi_profile = -R0 ./ v_safe .* grad_Vphi;
    
    % Reference mid-radius value for diagnostics
    [~, idx_mid] = min(abs(r_fine - a * 0.5));
    R_over_LVphi_mid = R_over_LVphi_profile(idx_mid);

    % =========================================================================
    % 2. GLOBAL PHYSICS PARAMETERS (MAGNITUDE)
    % =========================================================================
    % Ensure Global Maslov parameters exist (calculated in Stage 3)
    if isfield(derived_profiles, 'global_params') && isfield(derived_profiles.global_params, 'nu_star_e_solomon')
        nu_star_global = derived_profiles.global_params.nu_star_e_solomon;
        % Maslov's R/Ln (from linear fit 0.2-0.8) or Mid-radius fallback
        if isfield(derived_profiles.gradients, 'R_over_Ln_maslov')
            R_over_Ln_global = derived_profiles.gradients.R_over_Ln_maslov;
        else
            R_over_Ln_global = derived_profiles.gradients.R_over_Ln(idx_mid);
        end

        fprintf('   -> Global Nu*: %.4f\n', nu_star_global);
        fprintf('   -> Global R/Ln: %.4f\n', R_over_Ln_global);
        
        % Validity Check
        if nu_star_global > 0.8
            warning('TCABR:Physics', 'Global Nu* (%.2f) > 0.8. Solomon scaling is extrapolated.', nu_star_global);
        end
    else
        error('TCABR:Input', 'Global parameters missing. Run Stage 3 (Physics Profiles) first.');
    end

    % =========================================================================
    % 3. PHYSICS CONSTANTS & COEFFICIENTS
    % =========================================================================
    % Defined here to ensure scope availability for both MC loop and Table
    
    % Adjustable Physics Parameters
    F_gurcan = 1.0; % Asymmetry factor for TEP
    F_hahm   = 1.0; % Factor for Hahm (Usually 1.0, but kept explicit)

    % Empirical Coefficients (Solomon 2010)
    A_coef = 6.09;  A_err = 0.72;
    B_coef = 0.157; B_err = 0.072;
    C_coef = -24.2; C_err = 3.5;

    % =========================================================================
    % 4. MONTE CARLO SIMULATION
    % =========================================================================
    n_iter = constants.analysis.num_iterations;
    fprintf('... running %d MC iterations...\n', n_iter);
    
    % Sample Coefficients (Normal Distribution)
    A_s = A_coef + A_err * randn(n_iter, 1);
    B_s = B_coef + B_err * randn(n_iter, 1);
    C_s = C_coef + C_err * randn(n_iter, 1);

    % Storage
    chi_base_storage = zeros(n_iter, 1);
    models = {'Solomon', 'Hahm', 'Gurcan', 'Peeters_Rln2', 'Peeters_Rln_calc'};
    storage = struct();
    for m = 1:numel(models)
        storage.(models{m}).profile = zeros(n_iter, numel(r_fine));
    end

    % --- Calculation Loop ---
    for i = 1:n_iter
        % A. Base Diffusivity (Scalar)
        % chi = A*nu* + B*R/Ln. Uses Global parameters.
        chi_base = A_s(i) * nu_star_global + B_s(i) * R_over_Ln_global;
        chi_base_storage(i) = chi_base;

        % B. Pinch Velocities (Scalars)
        v_pinch.Solomon = C_s(i) * nu_star_global;
        v_pinch.Hahm    = -2 * chi_base / R0 * F_hahm;
        
        % Gürcan: -2(F + a/R0) + R/Ln (Global Geometry & Density)
        term_geo_gurc = -2 * (F_gurcan + epsilon_global);
        v_pinch.Gurcan  = (chi_base / R0) * (term_geo_gurc + R_over_Ln_global);
        
        % Peeters
        v_pinch.Peeters_Rln2     = (chi_base / R0) * (-4 - 2.0);
        v_pinch.Peeters_Rln_calc = (chi_base / R0) * (-4 - R_over_Ln_global);

        % C. Effective Diffusivity Profiles
        % chi_eff = chi_base * (1 + PinchNumber_Global / R_over_LVphi_Local(r))
        for m = 1:numel(models)
            name = models{m};
            vp = v_pinch.(name);
            PinchNumber = (R0 * vp) / chi_base;
            
            % Apply local geometry modulation
            storage.(name).profile(i,:) = chi_base .* (1 + PinchNumber ./ R_over_LVphi_profile');
        end
    end

    % =========================================================================
    % 5. POST-PROCESSING & PACKAGING
    % =========================================================================
    theory_results = struct();
    
    % A. Profiles (Mean & CI)
    for m = 1:numel(models)
        name = models{m};
        raw = storage.(name).profile;
        
        % Store directly under the model name (Variant: R0_profile)
        theory_results.(name).R0_profile.profile_avg = mean(raw, 1, 'omitnan')';
        ci = prctile(raw, [2.5, 97.5], 1);
        theory_results.(name).R0_profile.ci_lower = ci(1,:)';
        theory_results.(name).R0_profile.ci_upper = ci(2,:)';
    end

    % B. Metadata
    theory_results.meta.r_fine = r_fine;
    theory_results.meta.r_norm = r_norm;
    theory_results.meta.R_over_Ln_mid = R_over_Ln_global; % Used for plot labels
    
    % C. Base Diffusivity Statistics (for plotting the black dotted line)
    chi_base_mean = mean(chi_base_storage);
    chi_base_ci   = prctile(chi_base_storage, [2.5, 97.5]);
    
    theory_results.meta.chi_phi_solo_avg   = chi_base_mean * ones(size(r_fine));
    theory_results.meta.chi_phi_solo_lower = chi_base_ci(1) * ones(size(r_fine));
    theory_results.meta.chi_phi_solo_upper = chi_base_ci(2) * ones(size(r_fine));

    % =========================================================================
    % 6. DIAGNOSTIC TABLE (SCALAR VALUES)
    % =========================================================================
    % Calculate mean scalar values for the table using the mean coefficients
    chi_base_mean_tbl = A_coef * nu_star_global + B_coef * R_over_Ln_global;
    
    % Recompute mean pinch velocities for table display
    vp_solo_m = C_coef * nu_star_global;
    vp_hahm_m = -2 * chi_base_mean_tbl / R0 * F_hahm;
    vp_gurc_m = (chi_base_mean_tbl / R0) * (-2*(F_gurcan + epsilon_global) + R_over_Ln_global);
    vp_peet2_m = (chi_base_mean_tbl / R0) * (-4 - 2.0);
    vp_peetC_m = (chi_base_mean_tbl / R0) * (-4 - R_over_Ln_global);
    
    % Helper to calculate scalar effective diffusivity at r/a=0.5
    calc_eff = @(vp) chi_base_mean_tbl * (1 + (R0 * vp / chi_base_mean_tbl) / R_over_LVphi_mid);

    fprintf('\n----------------------------------------------------------------\n');
    fprintf('   GLOBAL TRANSPORT MODEL DIAGNOSTICS (Scalar Values)\n');
    fprintf('   Assumptions: <Te>=%.1f eV, <ne>=%.1e, R/Ln=%.2f, R/L_V=%.2f\n', ...
        derived_profiles.global_params.Te_avg, derived_profiles.global_params.ne_avg, ...
        R_over_Ln_global, R_over_LVphi_mid);
    fprintf('----------------------------------------------------------------\n');
    fprintf('| %-15s | %-10s | %-10s | %-10s |\n', 'Model', 'V_pinch', 'Chi_Base', 'Chi_Eff');
    fprintf('| %-15s | %-10s | %-10s | %-10s |\n', '', '[m/s]', '[m^2/s]', '[m^2/s]');
    fprintf('----------------------------------------------------------------\n');
    fprintf('| %-15s | %10.2f | %10.2f | %10.2f |\n', 'Solomon', vp_solo_m, chi_base_mean_tbl, calc_eff(vp_solo_m));
    fprintf('| %-15s | %10.2f | %10.2f | %10.2f |\n', 'Hahm (TEP)', vp_hahm_m, chi_base_mean_tbl, calc_eff(vp_hahm_m));
    fprintf('| %-15s | %10.2f | %10.2f | %10.2f |\n', 'Gurcan (TEP+n)', vp_gurc_m, chi_base_mean_tbl, calc_eff(vp_gurc_m));
    fprintf('| %-15s | %10.2f | %10.2f | %10.2f |\n', 'Peeters (R/Ln=2)', vp_peet2_m, chi_base_mean_tbl, calc_eff(vp_peet2_m));
    fprintf('| %-15s | %10.2f | %10.2f | %10.2f |\n', 'Peeters (Exp)', vp_peetC_m, chi_base_mean_tbl, calc_eff(vp_peetC_m));
    fprintf('----------------------------------------------------------------\n\n');

    fprintf('...theoretical variants calculated in %.2f s.\n', toc(t_start));
end