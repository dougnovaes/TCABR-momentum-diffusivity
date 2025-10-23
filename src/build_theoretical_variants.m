function [variants] = build_theoretical_variants(velocity_results, derived_profiles, constants, r_fine)
%BUILD_THEORETICAL_VARIANTS Systematically computes theoretical chi_phi_eff variants.
%   This function explores the main ambiguities in the formula for effective
%   momentum diffusivity from Peeters et al. (2007) by generating a matrix of
%   theoretical profiles. It combines different choices for the major radius
%   (global R0 vs. local R(r)) and the velocity gradient term (mid-radius
%   constant vs. full profile).
%
%   A Monte Carlo approach is used to propagate the uncertainties from the
%   empirical coefficients of the Solomon (2010) model, providing a mean
%   and 95% confidence interval for each variant.
%
%   Syntax:
%       variants = build_theoretical_variants(scaling_law_results, ...)

    arguments
        velocity_results (1,1) struct
        derived_profiles (1,1) struct
        constants (1,1) struct
        r_fine (:,1) {mustBeNumeric, mustBeReal, mustBeFinite}
    end

    fprintf('Building theoretical chi_eff variants with uncertainty propagation...\n');
    t_start = tic;

    % --- 1. Preparation and Robust Gradient Calculation ---
    a  = constants.machine.a;
    R0 = constants.machine.R0;
    
    % Define a safe radial mask for evaluation to avoid edge/core issues
    r_norm = r_fine / a;
    safe_mask = (r_norm >= 0.2) & (r_norm <= 0.9);

    % Robustly calculate the R/L_Vphi profile
    % Use a smoothed velocity profile to get a stable gradient
    v_phi_smooth = sgolayfilt(velocity_results.poly_fit_avg, 3, 15);
    grad_Vphi = gradient(v_phi_smooth, r_fine);
    
    % Regularize to avoid division by zero
    v_phi_safe = v_phi_smooth;
    v_phi_safe(abs(v_phi_safe) < 1e-3) = 1e-3; % Set a floor
    
    R_over_LVphi_profile = -R0 ./ v_phi_safe .* grad_Vphi;

    % Get mid-radius constant value
    [~, idx_mid] = min(abs(r_fine - a * 0.5));
    R_over_LVphi_mid = R_over_LVphi_profile(idx_mid);
    R_over_Ln_mid = derived_profiles.gradients.R_over_Ln(idx_mid);

    % --- 2. Monte Carlo Simulation for Uncertainty Propagation ---
    n_iter = 1000; % Number of MC iterations for CIs
    fprintf('...running %d Monte Carlo iterations for uncertainty bands.\n', n_iter);
    
    % Solomon et al. (2010) coefficients with errors
    A = 6.09;  A_err = 0.72;
    B = 0.157; B_err = 0.072;
    C = -24.2; C_err = 3.5;

    % Pre-sample coefficients
    A_samples = A + A_err * randn(n_iter, 1);
    B_samples = B + B_err * randn(n_iter, 1);
    C_samples = C + C_err * randn(n_iter, 1);
    
    nu_star_e = derived_profiles.collisionality.nu_star_e(:)'; % Ensure row vector for broadcasting
    R_coord = R0 * (1 + r_fine / R0);
    
    pinch_models = {'Solomon', 'Hahm', 'Gurcan', 'Peeters_Rln2', 'Peeters_Rln_calc'};
    
    % Pre-allocate storage for all profiles from all iterations
    storage = struct();
    for i = 1:numel(pinch_models)
        model_name = pinch_models{i};
        storage.(model_name).R0_mid         = zeros(n_iter, numel(r_fine));
        storage.(model_name).R0_profile     = zeros(n_iter, numel(r_fine));
        storage.(model_name).Rlocal_mid     = zeros(n_iter, numel(r_fine));
        storage.(model_name).Rlocal_profile = zeros(n_iter, numel(r_fine));
    end

    for i = 1:n_iter
        % Calculate base diffusivity for this iteration's sampled coefficients
        chi_phi_solo = A_samples(i) * nu_star_e + B_samples(i) * R_over_Ln_mid;
        
        % Calculate all pinch velocity profiles for this iteration
        v_pinch.Solomon = C_samples(i) * nu_star_e;
        v_pinch.Hahm = -2 * chi_phi_solo / R0;
        v_pinch.Gurcan = -(2 * chi_phi_solo ./ R_coord') .* (1 + r_fine' / R0);
        v_pinch.Peeters_Rln2 = (chi_phi_solo ./ R_coord') .* (-4 - 2);
        v_pinch.Peeters_Rln_calc = (chi_phi_solo ./ R_coord') .* (-4 - R_over_Ln_mid);
        
        % Calculate all 4 variants for each pinch model
        for m = 1:numel(pinch_models)
            model_name = pinch_models{m};
            vp = v_pinch.(model_name);
            
            % Variants using R0
            PinchNum_R0 = R0 * vp ./ chi_phi_solo;
            storage.(model_name).R0_mid(i,:)         = chi_phi_solo .* (1 + PinchNum_R0 / R_over_LVphi_mid);
            storage.(model_name).R0_profile(i,:)     = chi_phi_solo .* (1 + PinchNum_R0 ./ R_over_LVphi_profile');
            
            % Variants using R(r)
            PinchNum_Rlocal = R_coord' .* vp ./ chi_phi_solo;
            storage.(model_name).Rlocal_mid(i,:)     = chi_phi_solo .* (1 + PinchNum_Rlocal / R_over_LVphi_mid);
            storage.(model_name).Rlocal_profile(i,:) = chi_phi_solo .* (1 + PinchNum_Rlocal ./ R_over_LVphi_profile');
        end
    end
    
    % --- 3. Post-Process MC Results to get Mean and CIs ---
    variants = struct();
    for i = 1:numel(pinch_models)
        model_name = pinch_models{i};
        variant_names = fieldnames(storage.(model_name));
        for j = 1:numel(variant_names)
            variant_name = variant_names{j};
            
            profiles = storage.(model_name).(variant_name);
            
            % Calculate statistics and store in the final output structure
            variants.(model_name).(variant_name).profile_avg = mean(profiles, 1)';
            ci = prctile(profiles, [2.5, 97.5], 1);
            variants.(model_name).(variant_name).ci_lower = ci(1,:)';
            variants.(model_name).(variant_name).ci_upper = ci(2,:)';
        end
    end
    
    % --- 4. Package Metadata and Final Output ---
    variants.meta.r_fine = r_fine;
    variants.meta.r_norm = r_norm;
    variants.meta.safe_mask = safe_mask;
    variants.meta.R_over_LVphi_profile = R_over_LVphi_profile;
    variants.meta.R_over_LVphi_mid = R_over_LVphi_mid;
    variants.meta.R_over_Ln_mid = R_over_Ln_mid;
    variants.meta.chi_phi_solo_avg = mean(A_samples) * nu_star_e' + mean(B_samples) * R_over_Ln_mid;
    
    fprintf('...variant calculation finished in %.2f s.\n\n', toc(t_start));
end