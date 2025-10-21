function [diffusivity_results] = compute_diffusivity_profile(velocity_results, collision_profiles, r_fine)
%COMPUTE_DIFFUSIVITY_PROFILE Calculates the effective momentum diffusivity profile.
%   Produces profile_avg, ci_lower, ci_upper and compatibility aliases:
%   chi_eff_profile and chi_eff_error (half-width of CI).
%
%   Syntax:
%       diffusivity_results = compute_diffusivity_profile(velocity_results, collision_profiles, r_fine)

fprintf('Calculating effective momentum diffusivity profile...\n');

% --- 1. Extract necessary coefficients and profiles ---
chi_phi_coeff = velocity_results.chi_phi_coeff; % [m^2] scalar or vector
nu_iH0 = collision_profiles.nu_iH0;            % struct with .avg, .ci_lower, .ci_upper

% --- 2. Ensure shapes and types (column vectors, double) ---
r_fine = r_fine(:);
n_r = numel(r_fine);

% coefficient may be scalar or vector; force to double and conform
chi_phi_coeff = double(chi_phi_coeff);
if isscalar(chi_phi_coeff)
    chi_phi_coeff = chi_phi_coeff * ones(n_r,1);
else
    chi_phi_coeff = chi_phi_coeff(:);
    if numel(chi_phi_coeff) ~= n_r
        error('chi_phi_coeff length (%d) does not match r_fine (%d).', numel(chi_phi_coeff), n_r);
    end
end

% nu_iH0 fields: expect struct with .avg, .ci_lower, .ci_upper (vectors)
if ~isstruct(nu_iH0) || ~all(isfield(nu_iH0, {'avg','ci_lower','ci_upper'}))
    error('collision_profiles.nu_iH0 must be a struct with fields .avg, .ci_lower, .ci_upper');
end

nu_avg = double(nu_iH0.avg(:));
nu_lo  = double(nu_iH0.ci_lower(:));
nu_hi  = double(nu_iH0.ci_upper(:));

if any([numel(nu_avg), numel(nu_lo), numel(nu_hi)] ~= n_r)
    error('nu_iH0 vectors must have same length as r_fine (%d).', n_r);
end

% --- 3. Calculate the Diffusivity Profile and Propagate Uncertainty ---
chi_eff_avg      = chi_phi_coeff .* nu_avg;
chi_eff_ci_lower = chi_phi_coeff .* nu_lo;
chi_eff_ci_upper = chi_phi_coeff .* nu_hi;

% compute symmetric error (half-width) where meaningful
chi_eff_err = 0.5 .* (chi_eff_ci_upper - chi_eff_ci_lower);

% Guard against NaN/Inf from pathological inputs
chi_eff_avg(~isfinite(chi_eff_avg)) = 0;
chi_eff_ci_lower(~isfinite(chi_eff_ci_lower)) = 0;
chi_eff_ci_upper(~isfinite(chi_eff_ci_upper)) = 0;
chi_eff_err(~isfinite(chi_eff_err)) = 0;

% --- 4. Package Results into Output Structure (with compatibility aliases) ---
diffusivity_results.r_fine         = r_fine;
diffusivity_results.profile_avg    = chi_eff_avg;        % canonical name used in pipeline
diffusivity_results.ci_lower       = chi_eff_ci_lower;
diffusivity_results.ci_upper       = chi_eff_ci_upper;
diffusivity_results.chi_phi_coeff  = chi_phi_coeff;

% Compatibility: plots / utilities may expect these names
diffusivity_results.chi_eff_profile = chi_eff_avg;       % alias (old name)
diffusivity_results.chi_eff_error   = chi_eff_err;       % alias (half-width)
% also keep named CI fields for other code
diffusivity_results.ci_halfwidth    = chi_eff_err;       % extra alias

fprintf('Effective momentum diffusivity profile calculated successfully.\n\n');
end