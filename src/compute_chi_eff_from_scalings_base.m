function [scaling_law_results] = compute_chi_eff_from_scalings_base(velocity_results, derived_profiles, constants, r_fine)
%COMPUTE_CHI_EFF_FROM_SCALINGS_BASE Calculates effective diffusivity from theoretical scaling laws (Base Method).
%   This function computes theoretical profiles for chi_phi_eff based on the
%   formula from Peeters et al. (2007) (Eq. 3.14 in the thesis):
%       chi_phi_eff = chi_phi * (1 + (R*V_pinch/chi_phi) * (1 / (R/L_Vphi)))
%
%   This is the BASELINE implementation, consistent with the original analysis,
%   which uses characteristic, CONSTANT values for the normalised gradients
%   (R/L_n, R/L_Vphi) evaluated at mid-radius (r/a = 0.5).
%
%   Syntax:
%       scaling_law_results = compute_chi_eff_from_scalings_base(...)

    arguments
        velocity_results (1,1) struct
        derived_profiles (1,1) struct
        constants (1,1) struct
        r_fine (:,1) {mustBeNumeric, mustBeReal, mustBeFinite}
    end

    fprintf('Calculating effective diffusivity from scaling laws (Base Method)...\n');

    % --- 1. Extract Profiles and Parameters ---
    a         = constants.machine.a;
    R0        = constants.machine.R0;
    % **CORRECTION**: Use the specific Solomon-consistent collisionality
    nu_star_e = derived_profiles.collisionality.nu_star_e_solomon;

    % --- 2. Calculate Constant Gradient Values at Mid-Radius ---
    [~, idx_mid] = min(abs(r_fine - a * 0.5));

    % R/L_Vphi = -R0 * (1/V_phi) * (dV_phi/dr)
    Vphi_profile = velocity_results.poly_fit_avg;
    grad_Vphi = gradient(Vphi_profile, r_fine);
    R_over_LVphi_profile = -R0 ./ (Vphi_profile + eps) .* grad_Vphi;
    R_over_LVphi_const = R_over_LVphi_profile(idx_mid);

    R_over_Ln_const = derived_profiles.gradients.R_over_Ln(idx_mid);
    
    fprintf('... using characteristic values at r/a=0.5: R/L_n = %.2f, R/L_Vphi = %.2f\n', ...
        R_over_Ln_const, R_over_LVphi_const);

    % --- 3. Calculate Base Diffusivity and Pinch Velocity Profiles ---
    A = 6.09; B = 0.157; C = -24.2; % Solomon et al. 2010
    chi_phi_solo_profile = A * nu_star_e + B * R_over_Ln_const;

    v_pinch_solo   = C * nu_star_e;
    v_pinch_hahm   = -2 * chi_phi_solo_profile / R0;
    
    R_coord = R0 * (1 + r_fine / R0);
    F = 1.0; % Turbulence asymmetry factor
    v_pinch_gurcan = -(2 * chi_phi_solo_profile ./ R_coord) .* (F + r_fine / R0);
    
    v_pinch_peeters_Rln2   = (chi_phi_solo_profile ./ R_coord) .* (-4 - 2);
    v_pinch_peeters_Rln_calc = (chi_phi_solo_profile ./ R_coord) .* (-4 - R_over_Ln_const);

    % --- 4. Calculate Effective Diffusivity for Each Combination ---
    PinchNumber_Solo   = R0 * v_pinch_solo ./ chi_phi_solo_profile;
    PinchNumber_Hahm   = R0 * v_pinch_hahm ./ chi_phi_solo_profile;
    PinchNumber_Gurcan = R0 * v_pinch_gurcan ./ chi_phi_solo_profile;
    PinchNumber_Peet_2 = R0 * v_pinch_peeters_Rln2 ./ chi_phi_solo_profile;
    PinchNumber_Peet_C = R0 * v_pinch_peeters_Rln_calc ./ chi_phi_solo_profile;

    chi_eff_Solo   = chi_phi_solo_profile .* (1 + PinchNumber_Solo   / R_over_LVphi_const);
    chi_eff_Hahm   = chi_phi_solo_profile .* (1 + PinchNumber_Hahm   / R_over_LVphi_const);
    chi_eff_Gurcan = chi_phi_solo_profile .* (1 + PinchNumber_Gurcan / R_over_LVphi_const);
    chi_eff_Peet_2 = chi_phi_solo_profile .* (1 + PinchNumber_Peet_2 / R_over_LVphi_const);
    chi_eff_Peet_C = chi_phi_solo_profile .* (1 + PinchNumber_Peet_C / R_over_LVphi_const);

    % --- 5. Package Results ---
    scaling_law_results.r_fine = r_fine;
    scaling_law_results.chi_phi_solo = chi_phi_solo_profile;
    scaling_law_results.chi_eff.Solo = chi_eff_Solo;
    scaling_law_results.chi_eff.Hahm = chi_eff_Hahm;
    scaling_law_results.chi_eff.Gurcan = chi_eff_Gurcan;
    scaling_law_results.chi_eff.Peeters_Rln2 = chi_eff_Peet_2;
    scaling_law_results.chi_eff.Peeters_Rln_calc = chi_eff_Peet_C;
    scaling_law_results.meta.R_over_LVphi_const = R_over_LVphi_const;
    scaling_law_results.meta.R_over_Ln_const = R_over_Ln_const;

    fprintf('Calculation of base effective diffusivity from scaling laws complete.\n\n');
end