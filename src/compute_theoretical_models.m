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

fprintf('Calculating theoretical transport model profiles...\n');

% --- 1. Extract necessary variables ---
r_fine      = derived_profiles.r_fine;
a           = constants.machine.a;
B_pol       = magnetic_field.poloidal;
e_charge    = constants.physics.e_charge;
epsilon = constants.machine.epsilon_aspect_ratio;
nu_star_e   = derived_profiles.nu_star_e;
R_over_Ln   = derived_profiles.R_over_Ln;

% --- 2. Helander Model for Toroidal Velocity ---
% V_phi = (2*epsilon / (e*B_theta)) * (dT_i / dr)
% The dT_i/dr term is calculated analytically from the fitted canonical profile
% for maximum accuracy, avoiding numerical differentiation errors.
fprintf('... calculating Helander velocity profile\n');

% Extract the optimised parameters from the temperature fit
params_Ti = [temperature_results.params_optimized.T0; ...
             temperature_results.params_optimized.Ta; ...
             temperature_results.params_optimized.sigma];

% Analytical derivative of the canonical profile function T_i(r)
dTi_dr_analytical = analytical_dTi_dr(params_Ti, r_fine, a); % [eV/m]
delta = 0.0008; % Offset justified by Shafranov shift

% Use the physically-motivated offset for B_pol to avoid singularity
B_pol_offset  = B_pol + delta; % Offset justified by Shafranov shift

% Calculate the Helander velocity profile in m/s
% The factor of 'e_charge' is needed to convert d(Ti[eV])/dr to d(Ti[J])/dr
Vphi_Helander = (2 * epsilon ./ (e_charge * B_pol_offset)) .* (dTi_dr_analytical * e_charge); % [m/s]


% --- 3. Solomon Model for Diffusivity and Pinch ---
% Calculates momentum diffusivity (chi_phi) and pinch velocity (V_pinch)
% based on the empirical scaling laws from Solomon et al. (2010).
% chi_phi = A * nu_e* + B * R/Ln
% V_pinch = C * nu_e*
fprintf('... calculating Solomon diffusivity and pinch velocity profiles\n');

% Coefficients from Solomon et al. (2010), as cited in the thesis
A = 6.09; A_err = 0.72;
B = 0.157; B_err = 0.072;
C = -24.2; C_err = 3.5;

% --- Pinch Velocity (depends only on nu_star_e) ---
Vpinch_solo.profile = C * nu_star_e;
Vpinch_solo.upper_bound = (C + C_err) * nu_star_e;
Vpinch_solo.lower_bound = (C - C_err) * nu_star_e;

% --- Momentum Diffusivity (calculated for specific, constant R/Ln values) ---
% This follows the methodology in Solomon (2010) and the original thesis,
% where chi_phi is examined for characteristic values of R/Ln.

% Case 1: R/Ln = 2 (typical reference value)
chi_phi_solo_Rln2.profile = A * nu_star_e + B * 2;
chi_phi_solo_Rln2.upper_bound = (A + A_err) * nu_star_e + (B + B_err) * 2;
chi_phi_solo_Rln2.lower_bound = (A - A_err) * nu_star_e + (B - B_err) * 2;

% Case 2: R/Ln calculated at r/a = 0.5 (approx. 6 for TCABR)
% Find the index closest to r/a = 0.5
[~, idx_mid_radius] = min(abs(r_fine - a*0.5));
R_over_Ln_mid = R_over_Ln(idx_mid_radius);

% Correct profile and bounds (note the '+' for the B term)
chi_phi_solo_Rln_calc.profile = A * nu_star_e + B * R_over_Ln_mid;
chi_phi_solo_Rln_calc.upper_bound = (A + A_err) * nu_star_e + (B + B_err) * R_over_Ln_mid;
chi_phi_solo_Rln_calc.lower_bound = (A - A_err) * nu_star_e + (B - B_err) * R_over_Ln_mid;

% Ensure bounds are ordered (lower <= upper) element-wise
tmp_lower = min(chi_phi_solo_Rln_calc.lower_bound, chi_phi_solo_Rln_calc.upper_bound);
tmp_upper = max(chi_phi_solo_Rln_calc.lower_bound, chi_phi_solo_Rln_calc.upper_bound);
chi_phi_solo_Rln_calc.lower_bound = tmp_lower;
chi_phi_solo_Rln_calc.upper_bound = tmp_upper;

% --- 4. Package Results into Output Structure ---
theoretical_models.Vphi_Helander = Vphi_Helander;
theoretical_models.chi_phi_solo_Rln2 = chi_phi_solo_Rln2;
theoretical_models.chi_phi_solo_Rln_calc = chi_phi_solo_Rln_calc;
theoretical_models.R_over_Ln_mid_value = R_over_Ln_mid; % Store the value used
theoretical_models.Vpinch_solo = Vpinch_solo;

fprintf('Theoretical model profiles calculated successfully.\n\n');
end


% --- Local Helper Function for Analytical Derivative ---
function dTi_dr = analytical_dTi_dr(params, r, a)
    %ANALYTICAL_DTI_DR Calculates the analytical derivative of the canonical Ti profile.
    T0 = params(1);
    Ta = params(2);
    sigma = params(3);
    
    % Derivative of T_i(r) = (T0 - Ta) * (1 - (r/a)^2)^sigma + Ta
    base = 1 - (r ./ a).^2;
    base(base<0) = 0; % Prevent complex numbers
    
    dTi_dr = (T0 - Ta) .* sigma .* (base .^ (sigma - 1)) .* (-2 * r / a^2);
end