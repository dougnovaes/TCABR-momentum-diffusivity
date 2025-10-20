function [diffusivity_results] = compute_diffusivity_profile(velocity_results, collision_profiles, r_fine)
%COMPUTE_DIFFUSIVITY_PROFILE Calculates the effective momentum diffusivity profile.
%   This function computes the final effective momentum diffusivity, chi_phi_eff,
%   by combining the results from the velocity and collision analyses. It implements
%   the core relationship derived in the thesis (Eq. 3.13), where chi_phi_eff
%   is proportional to the ion-neutral collision frequency.
%
%   chi_phi_eff = [ (SUM(1/lambda_j^2)) * (3/(4*epsilon) - 1) ] * nu_iH0
%
%   The function also propagates the uncertainty from the nu_iH0 profile to the
%   final diffusivity profile.
%
%   Syntax:
%       diffusivity_results = compute_diffusivity_profile(velocity_results, collision_profiles, r_fine)
%
%   Inputs:
%       velocity_results   - Structure from analyze_velocity_profile, containing chi_phi_coeff.
%       collision_profiles - Structure from compute_collision_profiles, containing nu_iH0 profiles.
%       r_fine             - The high-resolution radial grid [m].
%
%   Output:
%       diffusivity_results - A structure containing the final mean chi_phi_eff profile
%                             and its associated confidence band.

fprintf('Calculating effective momentum diffusivity profile...\n');

% --- 1. Extract necessary coefficients and profiles ---
% The coefficient combining eigenvalues and geometry, calculated in the velocity analysis.
chi_phi_coeff = velocity_results.chi_phi_coeff; % [m^2]

% The ion-neutral collision frequency profiles, including uncertainty.
nu_iH0 = collision_profiles.nu_iH0; % [s^-1]

% --- 2. Calculate the Diffusivity Profile and Propagate Uncertainty ---
% Multiply the constant coefficient by the nu_iH0 profile and its confidence bands.
chi_eff_avg       = chi_phi_coeff * nu_iH0.avg;
chi_eff_ci_lower  = chi_phi_coeff * nu_iH0.ci_lower;
chi_eff_ci_upper  = chi_phi_coeff * nu_iH0.ci_upper;

% --- 3. Package Results into Output Structure ---
diffusivity_results.r_fine      = r_fine;
diffusivity_results.profile_avg = chi_eff_avg;
diffusivity_results.ci_lower    = chi_eff_ci_lower;
diffusivity_results.ci_upper    = chi_eff_ci_upper;
diffusivity_results.chi_phi_coeff = chi_phi_coeff; % Also store the coefficient for reference

fprintf('Effective momentum diffusivity profile calculated successfully.\n\n');

end