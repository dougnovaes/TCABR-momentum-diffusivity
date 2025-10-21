function nu_iH0 = compute_nu_iH0_from_nH0(n_H0, Ti_profiles_or_temp_results, ~, r_fine)
% COMPUTE_NU_IH0_FROM_NH0
%   Lightweight wrapper that produces ion-neutral collision frequency (nu_iH0)
%   from a neutral density profile n_H0(r). Internally it calls the existing
%   compute_collision_profiles(...) function so that all sanitisation and
%   percentile logic remains centralized.
%
% USAGE
%   nu_iH0 = compute_nu_iH0_from_nH0(n_H0, Ti_profiles_or_temp_results, constants, r_fine)
%
% INPUTS
%   n_H0                       - neutral density profile (n_r x 1) [m^-3]
%   Ti_profiles_or_temp_results- either:
%                                   * the full temperature_results struct (preferred), or
%                                   * a numeric Ti vector (n_r x 1) in eV, or
%                                   * a small struct with .profile_avg or .bootstrap.all_fitted_profiles
%   constants                  - project constants struct (passed to compute_collision_profiles)
%   r_fine                     - radial grid (n_r x 1) [m]
%
% OUTPUT
%   nu_iH0                     - struct with fields:
%                                   .avg      (n_r x 1)
%                                   .ci_lower (n_r x 1)
%                                   .ci_upper (n_r x 1)
%
% NOTES
%   - If only a single Ti profile is provided, the wrapper creates a
%     single-sample bootstrap input so that compute_collision_profiles
%     returns percentile fields (they will be equal to the mean).
%   - This wrapper intentionally does not duplicate the collision logic:
%     it adapts inputs and extracts nu_iH0 from compute_collision_profiles.
%
% Author: ChatGPT (adapted to user's codebase)
% Date: 2025

% ---- Basic checks and shape normalisation --------------------------------
if nargin < 4
    error('compute_nu_iH0_from_nH0 requires 4 inputs: n_H0, Ti_profiles_or_temp_results, constants, r_fine.');
end

n_H0 = double(n_H0(:));
r_fine = double(r_fine(:));
n_r = numel(r_fine);

if numel(n_H0) ~= n_r
    error('Length mismatch: n_H0 (length %d) must match r_fine (length %d).', numel(n_H0), n_r);
end

% ---- Prepare temperature_results in the expected form ----------------------
temp_in = Ti_profiles_or_temp_results;

% If user passed numeric Ti vector, create minimal temperature_results with one bootstrap sample
if isnumeric(temp_in) && isvector(temp_in)
    if numel(temp_in) ~= n_r
        error('Numeric Ti vector must have length equal to r_fine (%d).', n_r);
    end
    temperature_results = struct();
    % compute_collision_profiles expects bootstrap.all_fitted_profiles sized (iterations x points)
    temperature_results.bootstrap.all_fitted_profiles = double(temp_in(:)).'; % 1 x n_r
elseif isstruct(temp_in)
    temperature_results = temp_in;
    % If profile_avg present but no bootstrap, produce single sample bootstrap
    if isfield(temperature_results, 'profile_avg') && (~isfield(temperature_results,'bootstrap') || ~isfield(temperature_results.bootstrap,'all_fitted_profiles'))
        temperature_results.bootstrap.all_fitted_profiles = double(temperature_results.profile_avg(:)).';
    end
else
    error('Ti_profiles_or_temp_results must be a struct or numeric vector of length %d.', n_r);
end

% ---- Construct neutral_profile expected by compute_collision_profiles -----
neutral_profile = struct();
neutral_profile.n_H0 = n_H0(:);

% ---- Call compute_collision_profiles and extract nu_iH0 -------------------
try
    collision_profiles = compute_collision_profiles(temperature_results, neutral_profile, r_fine);
catch ME
    error('compute_collision_profiles failed. Message: %s', ME.message);
end

if ~isfield(collision_profiles, 'nu_iH0')
    error('compute_collision_profiles did not return collision_profiles.nu_iH0 as expected.');
end

nu_iH0 = collision_profiles.nu_iH0;

% Force column vectors and double type
nu_iH0.avg      = double(nu_iH0.avg(:));
nu_iH0.ci_lower = double(nu_iH0.ci_lower(:));
nu_iH0.ci_upper = double(nu_iH0.ci_upper(:));

% Final sanity check
if any([numel(nu_iH0.avg), numel(nu_iH0.ci_lower), numel(nu_iH0.ci_upper)] ~= n_r)
    error('Output nu_iH0 fields must each have length equal to r_fine (%d).', n_r);
end

end
