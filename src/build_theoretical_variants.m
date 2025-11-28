function [variants] = build_theoretical_variants(velocity_results, derived_profiles, constants, r_fine)
%BUILD_THEORETICAL_VARIANTS Systematically computes theoretical chi_phi_eff variants.
%   This function implements the "Hybrid Approach" for theoretical benchmarking:
%   1.  MAGNITUDE: Uses GLOBAL representative parameters (<Te>, <ne>, nu*_global)
%       to calculate the baseline turbulent diffusivity (chi_phi) and pinch velocity
%       (V_pinch). This avoids non-physical "explosions" at the cold edge.
%   2.  SHAPE: Uses the LOCAL experimental velocity gradient (R/L_Vphi) to 
%       calculate the effective diffusivity (chi_eff), preserving physical 
%       singularities where shear vanishes.
%
%   Syntax:
%       variants = build_theoretical_variants(velocity_results, derived_profiles, constants, r_fine)

    arguments
        velocity_results (1,1) struct
        derived_profiles (1,1) struct
        constants (1,1) struct
        r_fine (:,1) {mustBeNumeric, mustBeReal, mustBeFinite}
    end

    fprintf('Building theoretical chi_eff variants (Hybrid Global/Local)...\n');
    t_start = tic;

    % --- 1. Setup & Geometry ---
    a  = constants.machine.a;
    R0 = constants.machine.R0;
    r_norm = r_fine / a;
    R_coord = R0 * (1 + r_fine / R0); % Local Major Radius R(r)
    safe_mask = (r_norm >= 0.2) & (r_norm <= 0.9);

    % --- 2. Local Geometric Gradients (The "Shape") ---
    % Velocity Gradient (from polynomial fit)
    v_phi = velocity_results.poly_fit_avg; 
    grad_Vphi = gradient(v_phi, r_fine);
    v_safe = v_phi; v_safe(abs(v_safe) < 1e-9) = 1e-9; % Avoid div by zero crash
    
    % Profile of R/L_Vphi
    R_over_LVphi_profile = -R0 ./ v_safe .* grad_Vphi;

    % Mid-radius reference values
    [~, idx_mid] = min(abs(r_fine - a * 0.5));
    R_over_LVphi_mid = R_over_LVphi_profile(idx_mid);

    % R/L_n: Use fixed mid-radius value to stabilize Peeters/Solomon base
    if isfield(derived_profiles, 'R_over_Ln')
        R_over_Ln_mid = derived_profiles.R_over_Ln(idx_mid);
    else
        R_over_Ln_mid = derived_profiles.gradients.R_over_Ln(idx_mid);
    end
    
    % --- 3. Global Physics Parameters (The "Magnitude") ---
    % Retrieve the Global Solomon Collisionality calculated in Step 1
    if isfield(derived_profiles, 'global_params') && isfield(derived_profiles.global_params, 'nu_star_e_solomon')
        nu_star_global = derived_profiles.global_params.nu_star_e_solomon;
        fprintf('   Using Global Solomon Collisionality: %.4f\n', nu_star_global);
    else
        % Fallback if Step 1 wasn't re-run or structure differs
        warning('Global Solomon collisionality not found. Using mid-radius local value as fallback.');
        if isfield(derived_profiles, 'nu_star_e_solomon')
            nu_star_global = derived_profiles.nu_star_e_solomon(idx_mid);
        else
             nu_star_global = derived_profiles.nu_star_e(idx_mid);
        end
    end

    % --- 4. Constants & Monte Carlo Setup ---
    % Gürcan Asymmetry Factor (Exposed for tuning)
    F_gurcan = 0.5; % 0.5 = Mixed, 1.0 = Ballooning.

    % Solomon Coefficients (2010)
    A = 6.09;  A_err = 0.72;
    B = 0.157; B_err = 0.072;
    C = -24.2; C_err = 3.5;

    n_iter = 1000;
    fprintf('...running %d iterations...\n', n_iter);
    
    % Pre-sample coefficients
    A_s = A + A_err * randn(n_iter, 1);
    B_s = B + B_err * randn(n_iter, 1);
    C_s = C + C_err * randn(n_iter, 1);

    % Pre-allocate storage
    models = {'Solomon', 'Hahm', 'Gurcan', 'Peeters_Rln2', 'Peeters_Rln_calc'};
    storage = struct();
    for m = 1:numel(models)
        name = models{m};
        storage.(name).R0_profile = zeros(n_iter, numel(r_fine));  % Main variant (Global Base + Local Shape)
        storage.(name).R0_mid     = zeros(n_iter, numel(r_fine));  % Constant reference
        storage.(name).Rlocal_mid = zeros(n_iter, numel(r_fine));
        storage.(name).Rlocal_profile = zeros(n_iter, numel(r_fine));
    end

    % --- 5. Calculation Loop ---
    for i = 1:n_iter
        % A. Calculate SCALAR Base Diffusivity (Global Regime)
        % chi = A * nu_global + B * R/Ln_mid
        % This gives a stable "floor" for the diffusivity (~2-5 m^2/s)
        chi_base = A_s(i) * nu_star_global + B_s(i) * R_over_Ln_mid;

        % B. Calculate SCALAR Pinch Velocities (Global Regime)
        % Solomon: depends only on nu*
        vp_solo = C_s(i) * nu_star_global;

        % Prepare Pinch Velocity Struct (Scalar values, except Gürcan geometric factor)
        vp = struct();
        vp.Solomon = vp_solo;
        vp.Hahm    = -2 * chi_base / R0;
        
        % Gürcan: -(2*Chi/R) * (F + r/R0). 
        % The base magnitude is global, but the geometric factor (r/R0) is local.
        vp.Gurcan  = -(2 * chi_base ./ R_coord') .* (F_gurcan + r_fine'/R0);
        
        % Peeters: (Chi/R) * (-4 - R/Ln)
        % We use FIXED R/Ln to avoid edge noise in the pinch term magnitude
        vp.Peeters_Rln2     = (chi_base / R0) * (-6); % Fixed R/Ln=2 -> (-4-2) = -6
        vp.Peeters_Rln_calc = (chi_base / R0) * (-4 - R_over_Ln_mid); % Mid R/Ln

        % C. Calculate Effective Diffusivity Profiles
        % Formula: Chi_eff(r) = Chi_base * (1 + PinchNumber(r) / R_over_LVphi(r))
        
        for m = 1:numel(models)
            name = models{m};
            
            % Extract Pinch Velocity (Scalar or Vector)
            if strcmp(name, 'Gurcan')
                v_pinch = vp.Gurcan; % Vector
                PinchNum_R0 = (R0 * v_pinch) ./ chi_base;
                PinchNum_Rloc = (R_coord' .* v_pinch) ./ chi_base;
            else
                v_val = vp.(name); % Scalar
                PinchNum_R0 = (R0 * v_val) / chi_base;
                PinchNum_Rloc = (R_coord' * v_val) / chi_base;
            end

            % --- Variant 1: R0, Profile Grads (THE HYBRID MODEL) ---
            % Uses global chi/pinch but local Vphi gradients
            % This preserves the "singularities" where grad(V)=0 but stabilizes the magnitude.
            storage.(name).R0_profile(i,:) = chi_base .* (1 + PinchNum_R0 ./ R_over_LVphi_profile');

            % --- Variant 2: R0, Mid Grads (Constant Profile) ---
            % Reference line (flat)
            storage.(name).R0_mid(i,:) = chi_base .* (1 + PinchNum_R0 ./ R_over_LVphi_mid);

            % --- Variant 3/4: Local R (Minor variations) ---
            storage.(name).Rlocal_profile(i,:) = chi_base .* (1 + PinchNum_Rloc ./ R_over_LVphi_profile');
            storage.(name).Rlocal_mid(i,:) = chi_base .* (1 + PinchNum_Rloc ./ R_over_LVphi_mid);
        end
    end
    
    % --- 6. Package Statistics ---
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
    
    % Meta info
    variants.meta.r_fine = r_fine;
    variants.meta.r_norm = r_norm;
    variants.meta.safe_mask = safe_mask;
    variants.meta.R_over_LVphi_profile = R_over_LVphi_profile;
    variants.meta.R_over_LVphi_mid = R_over_LVphi_mid;
    variants.meta.R_over_Ln_mid = R_over_Ln_mid;
    
    % Save the mean global baseline diffusivity for plotting reference (constant line)
    chi_base_mean = mean(A_s * nu_star_global + B_s * R_over_Ln_mid);
    variants.meta.chi_phi_solo_avg = chi_base_mean * ones(size(r_fine));
    
    fprintf('...variants calculated in %.2f s.\n', toc(t_start));
end