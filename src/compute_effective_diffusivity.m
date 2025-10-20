function [eff_diffusivity] = compute_effective_diffusivity(velocity_results, collision_profiles, r_fine)
%COMPUTE_EFFECTIVE_DIFFUSIVITY Calculates the effective momentum diffusivity profile.
%   This function computes the effective momentum diffusivity, chi_phi_eff,
%   based on the method derived in the thesis (Eq. 3.13), which states:
%   chi_phi_eff = [ Sum(1/lambda_j^2) * (3/(4*epsilon) - 1) ] * nu_iH0
%
%   It combines the eigenvalue-derived coefficient from the velocity analysis
%   with the ion-neutral collision frequency profile. The function also
%   robustly propagates the uncertainty from the collision frequency profile
%   to the final diffusivity profile.
%
%   Syntax:
%       eff_diffusivity = compute_effective_diffusivity(velocity_results, collision_profiles, r_fine, constants)
%
%   Inputs:
%       velocity_results   - Structure from analyze_velocity_profile.
%       collision_profiles - Structure from compute_collision_profiles.
%       r_fine             - The high-resolution radial grid [m].
%       constants          - The main constants structure.
%
%   Output:
%       eff_diffusivity - A structure containing the mean profile and
%                         confidence band for chi_phi_eff.

fprintf('Calculating effective momentum diffusivity profile...\n');

% --- 1. Extract Necessary Data ---
% Coefficient derived from Bessel analysis: [ Sum(1/lambda_j^2) * (3/(4*epsilon) - 1) ]
% This term is a scalar constant.
chi_phi_coeff = velocity_results.chi_phi_coeff; % [m^2]

% Ion-neutral collision frequency profiles, including uncertainty bounds
nu_iH0_avg      = collision_profiles.nu_iH0.avg;
nu_iH0_ci_lower = collision_profiles.nu_iH0.ci_lower;
nu_iH0_ci_upper = collision_profiles.nu_iH0.ci_upper;

% --- 2. Calculate Diffusivity Profile and Propagate Uncertainty ---
% The uncertainty is propagated via simple multiplication of the constant
% coefficient with the mean and confidence bounds of the collision frequency.
profile_avg      = chi_phi_coeff * nu_iH0_avg;
profile_ci_lower = chi_phi_coeff * nu_iH0_ci_lower;
profile_ci_upper = chi_phi_coeff * nu_iH0_ci_upper;

% --- 3. Package Results into Output Structure ---
eff_diffusivity.r_fine = r_fine;
eff_diffusivity.profile_avg = profile_avg;
eff_diffusivity.ci_lower = profile_ci_lower;
eff_diffusivity.ci_upper = profile_ci_upper;

fprintf('Effective diffusivity profile calculated successfully.\n\n');

end