function [variants] = build_theoretical_variants(velocity_results, derived_profiles, constants, r_fine)
%BUILD_THEORETICAL_VARIANTS Computes theoretical transport variants using a Hybrid Approach.
%
%   PURPOSE:
%       Generates Monte Carlo uncertainty bands for theoretical momentum diffusivity (chi_phi)
%       and effective diffusivity (chi_eff) profiles.
%
%   METHODOLOGY (HYBRID GLOBAL/LOCAL):
%       This function bridges the gap between multi-machine scalings and local transport physics:
%       1. MAGNITUDE (Global): Uses the Global Solomon Collisionality (nu*_global) calculated
%          from volume-averaged parameters to set the baseline level of turbulent diffusivity.
%          This prevents unphysical divergences at the cold plasma edge.
%       2. SHAPE (Local): Uses the LOCAL experimental velocity gradient (R/L_Vphi) to 
%          modulate the effective diffusivity via the Pinch term. This preserves the 
%          physical topology of the transport, including singularities where shear vanishes.
%
%   INPUTS:
%       velocity_results : Struct with fitted V_phi and uncertainties.
%       derived_profiles : Struct with collisionality and density gradients.
%       constants        : Project constants.
%       r_fine           : Radial grid.
%
%   OUTPUTS:
%       variants         : Struct containing mean and CI profiles for:
%                          - Solomon (Base)
%                          - Hahm (TEP)
%                          - Gürcan (TEP + Density)
%                          - Peeters (Coriolis)

    arguments
        velocity_results (1,1) struct
        derived_profiles (1,1) struct
        constants (1,1) struct
        r_fine (:,1) {mustBeNumeric, mustBeReal, mustBeFinite}
    end

    fprintf('Building theoretical chi_eff variants (Hybrid Global/Local)...\n');
    t_start = tic;

    % =========================================================================
    % 1. SETUP & GEOMETRY
    % =========================================================================
    a  = constants.machine.a;
    R0 = constants.machine.R0;
    r_norm = r_fine / a;
    R_coord = R0 * (1 + r_fine / R0); 

    % =========================================================================
    % 2. LOCAL GRADIENTS (THE "SHAPE")
    % =========================================================================
    % Velocity Gradient (from polynomial fit)
    v_phi = velocity_results.poly_fit_avg; 
    grad_Vphi = gradient(v_phi, r_fine);
    v_safe = v_phi; 
    v_safe(abs(v_safe) < 1e-9) = 1e-9; 
    
    % Profile of R/L_Vphi
    R_over_LVphi_profile = -R0 ./ v_safe .* grad_Vphi;
    
    [~, idx_mid] = min(abs(r_fine - a * 0.5));
    R_over_LVphi_mid = R_over_LVphi_profile(idx_mid);

    % R/L_n: Use fixed mid-radius value for stability
    if isfield(derived_profiles, 'R_over_Ln')
        R_over_Ln_mid = derived_profiles.R_over_Ln(idx_mid);
        R_over_Ln_profile = derived_profiles.R_over_Ln;
    else
        R_over_Ln_mid = derived_profiles.gradients.R_over_Ln(idx_mid);
        R_over_Ln_profile = derived_profiles.gradients.R_over_Ln;
    end
    
    % =========================================================================
    % 3. GLOBAL PHYSICS PARAMETERS (THE "MAGNITUDE")
    % =========================================================================
    if isfield(derived_profiles, 'global_params') && isfield(derived_profiles.global_params, 'nu_star_e_solomon')
        nu_star_global = derived_profiles.global_params.nu_star_e_solomon;
        fprintf('   -> Using Global Solomon Collisionality: %.4f\n', nu_star_global);
        
        % Range Check
        if nu_star_global > 0.8
             warning('TCABR:Physics', ...
                ['Global Collisionality (%.2f) exceeds Solomon''s range (0.8). ' ...
                 'Theoretical diffusivity is an EXTRAPOLATION.'], nu_star_global);
        end
    else
        warning('Global collisionality not found. Using mid-radius fallback.');
        nu_star_global = derived_profiles.nu_star_e_solomon(idx_mid);
    end

    % =========================================================================
    % 4. MONTE CARLO SETUP
    % =========================================================================
    F_gurcan = 1.0; 
    F_hahm   = 1.0;

    % Solomon Coefficients (2010)
    A = 6.09;  A_err = 0.72;
    B = 0.157; B_err = 0.072;
    C = -24.2; C_err = 3.5;

    n_iter = 1000;
    fprintf('...running %d iterations...\n', n_iter);
    
    A_s = A + A_err * randn(n_iter, 1);
    B_s = B + B_err * randn(n_iter, 1);
    C_s = C + C_err * randn(n_iter, 1);

    chi_base_storage = zeros(n_iter, 1); 

    models = {'Solomon', 'Hahm', 'Gurcan', 'Peeters_Rln2', 'Peeters_Rln_calc'};
    storage = struct();
    for m = 1:numel(models)
        name = models{m};
        storage.(name).R0_profile = zeros(n_iter, numel(r_fine));
        storage.(name).R0_mid     = zeros(n_iter, numel(r_fine));
    end

    % =========================================================================
    % 5. CALCULATION LOOP
    % =========================================================================
    for i = 1:n_iter
        % A. Base Diffusivity (Scalar)
        chi_base = A_s(i) * nu_star_global + B_s(i) * R_over_Ln_mid;
        chi_base_storage(i) = chi_base;

        % B. Pinch Velocities
        vp_solo = C_s(i) * nu_star_global;
        
        vp = struct();
        vp.Solomon = vp_solo;
        vp.Hahm    = -2 * chi_base / R0 * F_hahm;
        
        % Gürcan: V = (chi/R) * [ -2(F + r/R) + R/Ln ]
        term_geometric = -2 * (F_gurcan + r_fine'/R0);
        term_density   = R_over_Ln_profile'; 
        vp.Gurcan      = (chi_base ./ R_coord') .* (term_geometric + term_density);
        
        % Peeters
        vp.Peeters_Rln2     = (chi_base / R0) * (-4 - 2);
        vp.Peeters_Rln_calc = (chi_base / R0) * (-4 - R_over_Ln_mid); 

        % C. Effective Diffusivity
        for m = 1:numel(models)
            name = models{m};
            
            if strcmp(name, 'Gurcan')
                v_pinch = vp.Gurcan; 
                PinchNum_R0 = (R0 * v_pinch) ./ chi_base;
            else
                v_val = vp.(name); 
                PinchNum_R0 = (R0 * v_val) / chi_base;
            end
            
            % Hybrid Variant
            storage.(name).R0_profile(i,:) = chi_base .* (1 + PinchNum_R0 ./ R_over_LVphi_profile');
            % Constant Reference Variant
            storage.(name).R0_mid(i,:)     = chi_base .* (1 + PinchNum_R0 ./ R_over_LVphi_mid);
        end
    end
    
    % =========================================================================
    % 6. PACKAGE RESULTS
    % =========================================================================
    variants = struct();
    for i = 1:numel(models)
        name = models{i};
        subfields = fieldnames(storage.(name));
        for j = 1:numel(subfields)
            sf = subfields{j};
            raw = storage.(name).(sf);
            variants.(name).(sf).profile_avg = mean(raw, 1, 'omitnan')';
            ci = prctile(raw, [2.5, 97.5], 1);
            variants.(name).(sf).ci_lower = ci(1,:)';
            variants.(name).(sf).ci_upper = ci(2,:)';
        end
    end
    
    % Metadata
    variants.meta.r_fine = r_fine;
    variants.meta.r_norm = r_norm;
    % Define mask here, close to where it's stored
    variants.meta.safe_mask = (r_norm >= 0.2) & (r_norm <= 0.9); 
    variants.meta.R_over_LVphi_profile = R_over_LVphi_profile;
    variants.meta.R_over_LVphi_mid = R_over_LVphi_mid;
    variants.meta.R_over_Ln_mid = R_over_Ln_mid;
    
    % Base Statistics
    chi_base_mean = mean(chi_base_storage);
    chi_base_ci   = prctile(chi_base_storage, [2.5, 97.5]);
    
    variants.meta.chi_phi_solo_avg   = chi_base_mean * ones(size(r_fine));
    variants.meta.chi_phi_solo_lower = chi_base_ci(1) * ones(size(r_fine));
    variants.meta.chi_phi_solo_upper = chi_base_ci(2) * ones(size(r_fine));
    
    fprintf('...variants calculated in %.2f s.\n', toc(t_start));
end