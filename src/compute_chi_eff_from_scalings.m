function [scaling_law_results] = compute_chi_eff_from_scalings(velocity_results, derived_profiles, constants, r_fine)
%COMPUTE_CHI_EFF_FROM_SCALINGS Calculates effective diffusivity from theoretical scaling laws.
%   This function computes theoretical profiles for the effective momentum
%   diffusivity (chi_phi_eff) based on the formula from Peeters et al. (2007)
%   (Eq. 3.14 in the thesis):
%       chi_phi_eff = chi_phi * (1 + (R*V_pinch/chi_phi) * (1 / (R/L_Vphi)))
%
%   It uses a single model for the 'pure' diffusivity (Solomon) and combines
%   it with various theoretical models for the pinch velocity (Solomon, Hahm,
%   Gürcan, Peeters) to generate a set of comparable chi_phi_eff profiles.
%
%   Syntax:
%       scaling_law_results = compute_chi_eff_from_scalings(velocity_results, derived_profiles, constants, r_fine)
%
%   Inputs:
%       velocity_results - Structure from analyze_velocity_profile.
%       derived_profiles - Structure with nu_star_e, R_over_Ln, etc.
%       constants        - Structure with all project constants.
%       r_fine           - The high-resolution radial grid [m].
%
%   Output:
%       scaling_law_results - A structure containing all calculated theoretical
%                             profiles and intermediate values.

fprintf('Calculating effective diffusivity from scaling laws...\n');

% --- 1. Extract Necessary Profiles and Parameters ---
a         = constants.machine.a;
R0        = constants.machine.R0;
nu_star_e = derived_profiles.nu_star_e;

% --- 2. Calculate Constant Values for R/L_n and R/L_Vphi ---
% As per the methodology, we use characteristic values for the normalised
% gradients evaluated at mid-radius (r/a = 0.5) instead of the full profiles.
[~, idx_mid] = min(abs(r_fine - a * 0.5));

% Calculate R/L_Vphi = -R0 * (1/V_phi) * (dV_phi/dr)
Vphi_profile = velocity_results.poly_fit_avg;
grad_Vphi = gradient(Vphi_profile, r_fine);
R_over_LVphi_profile = -R0 ./ Vphi_profile .* grad_Vphi;
R_over_LVphi_const = R_over_LVphi_profile(idx_mid); % Use value at r/a=0.5

% Get the value of R/L_n at r/a = 0.5 for the Solomon diffusivity calculation
R_over_Ln_const_calc = derived_profiles.R_over_Ln(idx_mid);

fprintf('... using characteristic values at r/a=0.5: R/L_n = %.2f, R/L_Vphi = %.2f\n', ...
    R_over_Ln_const_calc, R_over_LVphi_const);

% --- 3. Calculate Base Diffusivity and Pinch Velocity Profiles ---
% These are the building blocks for the final chi_phi_eff calculations.

% Solomon momentum diffusivity, chi_phi^(Solo) (using constant R/L_n)
% Coefficients from Solomon et al. (2010)
A = 6.09; B = 0.157; C = -24.2;
chi_phi_solo_profile = A * nu_star_e + B * R_over_Ln_const_calc;

% Solomon pinch velocity, V_pinch^(Solo)
v_pinch_solo_profile = C * nu_star_e;

% Hahm pinch velocity, V_pinch^(Hahm) = -2*chi_phi/R0
v_pinch_hahm_profile = -2 * chi_phi_solo_profile / R0;

% Gürcan pinch velocity, V_pinch^(Gürc) = -(2*chi_phi/R)*(F + r/R0)
% For now, F (turbulence asymmetry factor) is assumed to be 1.
F = 1.0;
R_coord = R0 * (1 + r_fine / R0);
v_pinch_gurcan_profile = -(2 * chi_phi_solo_profile ./ R_coord) .* (F + r_fine / R0);

% Peeters pinch velocity, V_pinch^(Peet) = (chi_phi/R)*(-4 - R/L_n)
% We calculate this for two cases as done in the thesis.
v_pinch_peeters_Rln2_profile   = (chi_phi_solo_profile ./ R_coord) .* (-4 - 2);
v_pinch_peeters_Rln_calc_profile = (chi_phi_solo_profile ./ R_coord) .* (-4 - R_over_Ln_const_calc);


% --- 4. Calculate Effective Diffusivity for Each Combination ---
% Apply Eq. 3.14: chi_eff = chi_phi * (1 + PinchNumber / (R/L_Vphi))
fprintf('... combining models to compute effective diffusivity profiles\n');

PinchNumber_Solo   = R0 * v_pinch_solo_profile ./ chi_phi_solo_profile;
PinchNumber_Hahm   = R0 * v_pinch_hahm_profile ./ chi_phi_solo_profile;
PinchNumber_Gurcan = R0 * v_pinch_gurcan_profile ./ chi_phi_solo_profile;
PinchNumber_Peet_2 = R0 * v_pinch_peeters_Rln2_profile ./ chi_phi_solo_profile;
PinchNumber_Peet_C = R0 * v_pinch_peeters_Rln_calc_profile ./ chi_phi_solo_profile;

chi_eff_Solo   = chi_phi_solo_profile .* (1 + PinchNumber_Solo   / R_over_LVphi_const);
chi_eff_Hahm   = chi_phi_solo_profile .* (1 + PinchNumber_Hahm   / R_over_LVphi_const);
chi_eff_Gurcan = chi_phi_solo_profile .* (1 + PinchNumber_Gurcan / R_over_LVphi_const);
chi_eff_Peet_2 = chi_phi_solo_profile .* (1 + PinchNumber_Peet_2 / R_over_LVphi_const);
chi_eff_Peet_C = chi_phi_solo_profile .* (1 + PinchNumber_Peet_C / R_over_LVphi_const);

% --- 5. Package Results into Output Structure ---
scaling_law_results.r_fine = r_fine;
scaling_law_results.chi_phi_solo = chi_phi_solo_profile;
scaling_law_results.v_pinch_solo = v_pinch_solo_profile;
scaling_law_results.v_pinch_hahm = v_pinch_hahm_profile;
scaling_law_results.v_pinch_gurcan = v_pinch_gurcan_profile;
scaling_law_results.v_pinch_peeters_Rln2 = v_pinch_peeters_Rln2_profile;
scaling_law_results.v_pinch_peeters_Rln_calc = v_pinch_peeters_Rln_calc_profile;
scaling_law_results.chi_eff_Solo = chi_eff_Solo;
scaling_law_results.chi_eff_Hahm = chi_eff_Hahm;
scaling_law_results.chi_eff_Gurcan = chi_eff_Gurcan;
scaling_law_results.chi_eff_Peeters_Rln2 = chi_eff_Peet_2;
scaling_law_results.chi_eff_Peeters_Rln_calc = chi_eff_Peet_C;
scaling_law_results.R_over_LVphi_const = R_over_LVphi_const;
scaling_law_results.R_over_Ln_const = R_over_Ln_const_calc;

fprintf('Calculation of effective diffusivity from scaling laws complete.\n\n');

end