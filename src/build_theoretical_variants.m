function [variants] = build_theoretical_variants(velocity_results, derived_profiles, constants, r_fine)
%BUILD_THEORETICAL_VARIANTS Systematically computes theoretical chi_phi_eff variants.
%   This function explores the main ambiguities in the formula for effective
%   momentum diffusivity from Peeters et al. (2007) by generating a matrix of
%   theoretical profiles.
%
%   A Monte Carlo approach is used to propagate the uncertainties from the
%   empirical coefficients of the Solomon (2010) model.
%
%   NOTE: This version performs DIRECT CALCULATION using Solomon-consistent 
%   collisionality and local gradients. Singularities where R/L_Vphi -> 0 
%   are expected.

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
    
    r_norm = r_fine / a;
    safe_mask = (r_norm >= 0.2) & (r_norm <= 0.9);

    % Calculate Gradient directly from the polynomial fit (already smooth)
    v_phi_profile = velocity_results.poly_fit_avg; 
    grad_Vphi = gradient(v_phi_profile, r_fine);
    
    % Simple guard for absolute zero in Vphi to avoid immediate crash
    v_phi_safe = v_phi_profile;
    v_phi_safe(abs(v_phi_safe) < 1e-9) = 1e-9;
    
    R_over_LVphi_profile = -R0 ./ v_phi_safe .* grad_Vphi;

    % Get mid-radius constant value
    [~, idx_mid] = min(abs(r_fine - a * 0.5));
    R_over_LVphi_mid = R_over_LVphi_profile(idx_mid);
    
    % **CORRECTION 1**: Access R_over_Ln directly (flat structure)
    if isfield(derived_profiles, 'R_over_Ln')
        R_over_Ln_mid = derived_profiles.R_over_Ln(idx_mid);
    elseif isfield(derived_profiles, 'gradients') && isfield(derived_profiles.gradients, 'R_over_Ln')
        R_over_Ln_mid = derived_profiles.gradients.R_over_Ln(idx_mid);
    else
        error('R_over_Ln field not found in derived_profiles.');
    end

    % --- 2. Monte Carlo Simulation ---
    n_iter = 1000;
    fprintf('...running %d Monte Carlo iterations for uncertainty bands.\n', n_iter);
    
    % Coefficients (Solomon 2010)
    A = 6.09;  A_err = 0.72;
    B = 0.157; B_err = 0.072;
    C = -24.2; C_err = 3.5;

    A_samples = A + A_err * randn(n_iter, 1);
    B_samples = B + B_err * randn(n_iter, 1);
    C_samples = C + C_err * randn(n_iter, 1);
    
    % **CORRECTION 2**: Use Solomon-consistent collisionality
    if isfield(derived_profiles, 'nu_star_e_solomon')
        nu_star_e = derived_profiles.nu_star_e_solomon(:)';
    elseif isfield(derived_profiles.collisionality, 'nu_star_e_solomon')
         nu_star_e = derived_profiles.collisionality.nu_star_e_solomon(:)';
    else
        warning('Solomon collisionality not found in build_variants. Using Wesson fallback.');
        nu_star_e = derived_profiles.nu_star_e(:)';
    end
    
    R_coord = R0 * (1 + r_fine / R0);
    pinch_models = {'Solomon', 'Hahm', 'Gurcan', 'Peeters_Rln2', 'Peeters_Rln_calc'};
    
    % Pre-allocate
    storage = struct();
    for i = 1:numel(pinch_models)
        model_name = pinch_models{i};
        storage.(model_name).R0_mid         = zeros(n_iter, numel(r_fine));
        storage.(model_name).R0_profile     = zeros(n_iter, numel(r_fine));
        storage.(model_name).Rlocal_mid     = zeros(n_iter, numel(r_fine));
        storage.(model_name).Rlocal_profile = zeros(n_iter, numel(r_fine));
    end

    for i = 1:n_iter
        % Base diffusivity
        chi_phi_solo = A_samples(i) * nu_star_e + B_samples(i) * R_over_Ln_mid;
        
        % Pinch velocities
        v_pinch.Solomon = C_samples(i) * nu_star_e;
        v_pinch.Hahm = -2 * chi_phi_solo / R0;
        v_pinch.Gurcan = -(2 * chi_phi_solo ./ R_coord') .* (0.5 + r_fine' / R0);
        v_pinch.Peeters_Rln2 = (chi_phi_solo ./ R_coord') .* (-4 - 2);
        v_pinch.Peeters_Rln_calc = (chi_phi_solo ./ R_coord') .* (-4 - R_over_Ln_mid);
        
        for m = 1:numel(pinch_models)
            model_name = pinch_models{m};
            vp = v_pinch.(model_name);
            
            % 1. Global R0
            PinchNum_R0 = R0 * vp ./ chi_phi_solo;
            storage.(model_name).R0_mid(i,:)     = chi_phi_solo .* (1 + PinchNum_R0 / R_over_LVphi_mid);
            storage.(model_name).R0_profile(i,:) = chi_phi_solo .* (1 + PinchNum_R0 ./ R_over_LVphi_profile');
            
            % 2. Local R(r)
            PinchNum_Rlocal = R_coord' .* vp ./ chi_phi_solo;
            storage.(model_name).Rlocal_mid(i,:)     = chi_phi_solo .* (1 + PinchNum_Rlocal / R_over_LVphi_mid);
            storage.(model_name).Rlocal_profile(i,:) = chi_phi_solo .* (1 + PinchNum_Rlocal ./ R_over_LVphi_profile');
        end
    end
    
    % --- 3. Post-Process ---
    variants = struct();
    for i = 1:numel(pinch_models)
        model_name = pinch_models{i};
        variant_names = fieldnames(storage.(model_name));
        for j = 1:numel(variant_names)
            variant_name = variant_names{j};
            profiles = storage.(model_name).(variant_name);
            
            variants.(model_name).(variant_name).profile_avg = mean(profiles, 1, 'omitnan')';
            ci = prctile(profiles, [2.5, 97.5], 1);
            variants.(model_name).(variant_name).ci_lower = ci(1,:)';
            variants.(model_name).(variant_name).ci_upper = ci(2,:)';
        end
    end
    
    % --- 4. Package ---
    variants.meta.r_fine = r_fine;
    variants.meta.r_norm = r_norm;
    variants.meta.safe_mask = safe_mask;
    variants.meta.R_over_LVphi_profile = R_over_LVphi_profile;
    variants.meta.R_over_LVphi_mid = R_over_LVphi_mid;
    variants.meta.R_over_Ln_mid = R_over_Ln_mid;
    variants.meta.chi_phi_solo_avg = mean(A_samples) * nu_star_e' + mean(B_samples) * R_over_Ln_mid;
    
    fprintf('...variant calculation finished in %.2f s.\n\n', toc(t_start));
end