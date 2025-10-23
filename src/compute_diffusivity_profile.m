function [diffusivity_results] = compute_diffusivity_profile(velocity_results, collision_profiles, r_fine)
%COMPUTE_DIFFUSIVITY_PROFILE Calculates the effective momentum diffusivity profile ("Thesis Method").
%   This function computes the effective momentum diffusivity, chi_phi_eff,
%   based on the method derived in the thesis (Eq. 3.13):
%   chi_phi_eff = [ Sum(1/lambda_j^2) * (3/(4*epsilon) - 1) ] * nu_iH0
%
%   It combines the scalar, eigenvalue-derived coefficient from the velocity analysis
%   with the ion-neutral collision frequency profile, robustly propagating
%   the uncertainty from the collision frequency.
%
%   Syntax:
%       diffusivity_results = compute_diffusivity_profile(velocity_results, collision_profiles, r_fine)

arguments
        velocity_results (1,1) struct
        collision_profiles (1,1) struct
        r_fine (:,1) {mustBeNumeric, mustBeReal, mustBeFinite}
    end

    fprintf('Calculating effective momentum diffusivity profile (Thesis Method)...\n');

    % --- 1. Input Validation and Data Extraction ---
    if ~isfield(velocity_results, 'chi_phi_coeff')
        error('Input struct ''velocity_results'' is missing required field: ''chi_phi_coeff''');
    end
    if ~isfield(collision_profiles, 'nu_iH0')
        error('Input struct ''collision_profiles'' is missing required field: ''nu_iH0''');
    end

    chi_phi_coeff = velocity_results.chi_phi_coeff;
    nu_iH0 = collision_profiles.nu_iH0;
    
    validateattributes(chi_phi_coeff, {'numeric'}, {'scalar', 'real', 'finite'}, mfilename, 'chi_phi_coeff');
    if ~isfield(nu_iH0, 'avg') || ~isfield(nu_iH0, 'ci_lower') || ~isfield(nu_iH0, 'ci_upper')
         error('Input struct ''collision_profiles.nu_iH0'' is missing required fields.');
    end

    % --- 2. Calculate Diffusivity Profile and Propagate Uncertainty ---
    chi_eff_avg      = chi_phi_coeff * nu_iH0.avg;
    chi_eff_ci_lower = chi_phi_coeff * nu_iH0.ci_lower;
    chi_eff_ci_upper = chi_phi_coeff * nu_iH0.ci_upper;
    chi_eff_err      = 0.5 * (chi_eff_ci_upper - chi_eff_ci_lower);

    % --- 3. Package Results ---
    diffusivity_results.r_fine      = r_fine;
    diffusivity_results.profile_avg = chi_eff_avg;
    diffusivity_results.ci_lower    = chi_eff_ci_lower;
    diffusivity_results.ci_upper    = chi_eff_ci_upper;
    diffusivity_results.ci_halfwidth = chi_eff_err;
    diffusivity_results.chi_phi_coeff = chi_phi_coeff;

    % --- Compatibility Aliases (for older plotting scripts if needed) ---
    % Note: The canonical names are .profile_avg, .ci_lower, .ci_upper
    diffusivity_results.chi_eff_profile = chi_eff_avg;
    diffusivity_results.chi_eff_error   = chi_eff_err;

    fprintf('Effective diffusivity profile calculated successfully.\n\n');
end