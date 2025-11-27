function [theoretical_models] = compute_theoretical_models(temperature_results, magnetic_field, derived_profiles, constants)
%COMPUTE_THEORETICAL_MODELS Calculates profiles from theoretical transport models.
%   This function uses previously calculated plasma profiles to compute
%   theoretical predictions for key transport quantities. It currently
%   implements the Helander model for toroidal velocity and the Solomon
%   model for momentum diffusivity and pinch velocity.
%
%   Syntax:
%       theoretical_models = compute_theoretical_models(temp_results, mag_field, derived_profiles, constants)
%
%   Inputs:
%       temperature_results - Structure with fitted Ti profile and parameters.
%       magnetic_field     - Structure with B-field profiles.
%       derived_profiles   - Structure with nu_star, R/Ln, etc.
%       constants          - Structure with all project constants.
%
%   Output:
%       theoretical_models - A structure containing all calculated theoretical profiles.

    arguments
        temperature_results (1,1) struct
        magnetic_field (1,1) struct
        derived_profiles (1,1) struct
        constants (1,1) struct
    end

    fprintf('Calculating theoretical transport model profiles...\n');

    % --- 1. Extract necessary variables ---
    r_fine      = derived_profiles.r_fine;
    a           = constants.machine.a;
    B_pol       = magnetic_field.poloidal;
    e_charge    = constants.physics.e_charge;
    epsilon     = constants.machine.epsilon_aspect_ratio;
    
    % USE SOLOMON COLLISIONALITY (Drift based) for correct physics comparison
    if isfield(derived_profiles, 'nu_star_e_solomon')
        nu_star_e = derived_profiles.nu_star_e_solomon;
    elseif isfield(derived_profiles.collisionality, 'nu_star_e_solomon')
        nu_star_e = derived_profiles.collisionality.nu_star_e_solomon;
    else
        warning('Solomon collisionality not found. Using standard Wesson collisionality as fallback.');
        nu_star_e = derived_profiles.nu_star_e;
    end

    % Access R_over_Ln directly (flat structure)
    if isfield(derived_profiles, 'R_over_Ln')
        R_over_Ln = derived_profiles.R_over_Ln;
    else
        error('R_over_Ln field not found in derived_profiles structure.');
    end

    % --- 2. Helander Model for Toroidal Velocity ---
    % V_phi = (2*epsilon / (e*B_theta)) * (dT_i / dr)
    fprintf('... calculating Helander velocity profile\n');

    % Extract the optimised parameters from the temperature fit
    params_Ti = [temperature_results.params_optimized.T0; ...
                 temperature_results.params_optimized.Ta; ...
                 temperature_results.params_optimized.sigma];

    % Analytical derivative of the canonical profile function T_i(r)
    dTi_dr_analytical = analytical_dTi_dr(params_Ti, r_fine, a); % [eV/m]
    delta = 0.0008; % Offset justified by Shafranov shift

    % Use the physically-motivated offset for B_pol to avoid singularity
    B_pol_offset  = B_pol + delta; 

    % Calculate the Helander velocity profile in m/s
    Vphi_Helander = (2 * epsilon ./ (e_charge * B_pol_offset)) .* (dTi_dr_analytical * e_charge); 


    % --- 3. Solomon Model for Diffusivity and Pinch ---
    % Calculates momentum diffusivity (chi_phi) and pinch velocity (V_pinch)
    % based on the empirical scaling laws from Solomon et al. (2010).
    fprintf('... calculating Solomon diffusivity and pinch velocity profiles\n');

    % Coefficients from Solomon et al. (2010)
    A = 6.09; A_err = 0.72;
    B = 0.157; B_err = 0.072;
    C = -24.2; C_err = 3.5;

    % --- Pinch Velocity (depends only on nu_star_e) ---
    Vpinch_solo.profile = C * nu_star_e;
    Vpinch_solo.upper_bound = (C + C_err) * nu_star_e;
    Vpinch_solo.lower_bound = (C - C_err) * nu_star_e;

    % --- Momentum Diffusivity ---
    
    % Case 1: R/Ln = 2 (typical reference value)
    chi_phi_solo_Rln2.profile = A * nu_star_e + B * 2;
    chi_phi_solo_Rln2.upper_bound = (A + A_err) * nu_star_e + (B + B_err) * 2;
    chi_phi_solo_Rln2.lower_bound = (A - A_err) * nu_star_e + (B - B_err) * 2;

    % Case 2: R/Ln calculated at r/a = 0.5
    [~, idx_mid] = min(abs(r_fine - a*0.5));
    R_over_Ln_mid = R_over_Ln(idx_mid);

    chi_phi_solo_Rln_calc.profile = A * nu_star_e + B * R_over_Ln_mid;
    chi_phi_solo_Rln_calc.upper_bound = (A + A_err) * nu_star_e + (B + B_err) * R_over_Ln_mid;
    chi_phi_solo_Rln_calc.lower_bound = (A - A_err) * nu_star_e + (B - B_err) * R_over_Ln_mid;

    % Ensure bounds are ordered
    tmp_lower = min(chi_phi_solo_Rln_calc.lower_bound, chi_phi_solo_Rln_calc.upper_bound);
    tmp_upper = max(chi_phi_solo_Rln_calc.lower_bound, chi_phi_solo_Rln_calc.upper_bound);
    chi_phi_solo_Rln_calc.lower_bound = tmp_lower;
    chi_phi_solo_Rln_calc.upper_bound = tmp_upper;

    % --- 4. Package Results ---
    theoretical_models.Vphi_Helander = Vphi_Helander;
    theoretical_models.chi_phi_solo_Rln2 = chi_phi_solo_Rln2;
    theoretical_models.chi_phi_solo_Rln_calc = chi_phi_solo_Rln_calc;
    theoretical_models.R_over_Ln_mid_value = R_over_Ln_mid;
    theoretical_models.Vpinch_solo = Vpinch_solo;

    fprintf('Theoretical model profiles calculated successfully.\n\n');
end

function dTi_dr = analytical_dTi_dr(params, r, a)
    %ANALYTICAL_DTI_DR Calculates the analytical derivative of the canonical Ti profile.
    T0 = params(1); Ta = params(2); sigma = params(3);
    base = 1 - (r ./ a).^2;
    base(base<0) = 0;
    dTi_dr = (T0 - Ta) .* sigma .* (base .^ (sigma - 1)) .* (-2 * r / a^2);
end