function [collision_profiles] = compute_collision_profiles(temperature_results, neutral_profile, r_fine)
%COMPUTE_COLLISION_PROFILES Calculates charge exchange rate and ion-neutral collision frequency robustly.
%   This function computes the charge exchange (cx) rate coefficient and the
%   ion-neutral collision frequency (nu_iH0). It is designed to be highly
%   robust, with extensive input validation, error handling, and defensive
%   checks to ensure numerical stability and physical validity.
%
%   Syntax:
%       collision_profiles = compute_collision_profiles(temperature_results, neutral_profile, r_fine)

    arguments
        temperature_results (1,1) struct
        neutral_profile (1,1) struct
        r_fine (:,1) {mustBeNumeric, mustBeReal, mustBeFinite}
    end

    fprintf('Calculating collision rate and frequency profiles...\n');

    % --- 1. Input Validation and Data Sanitisation ---
    % Validate temperature_results structure
    if ~isfield(temperature_results, 'bootstrap') || ~isfield(temperature_results.bootstrap, 'all_fitted_profiles')
        error('Input structure ''temperature_results'' is missing the field: bootstrap.all_fitted_profiles');
    end
    
    % Validate neutral_profile structure
    if ~isfield(neutral_profile, 'n_H0')
        error('Input structure ''neutral_profile'' is missing the field: n_H0');
    end

    all_ti_profiles = temperature_results.bootstrap.all_fitted_profiles;

    % Ensure the matrix is 2D and has the expected orientation (iterations x points)
    if ~ismatrix(all_ti_profiles)
        error('Input ''all_fitted_profiles'' must be a 2D matrix.');
    end
    
    % Handle orientation if necessary
    if size(all_ti_profiles, 2) ~= numel(r_fine)
        if size(all_ti_profiles, 1) == numel(r_fine)
            all_ti_profiles = all_ti_profiles'; % Transpose to (iterations x points)
        else
            error('Dimension mismatch: The number of points in temperature profiles does not match the radial grid.');
        end
    end

    % Enforce physical and numerical validity (Clamping)
    all_ti_profiles = real(double(all_ti_profiles)); % Discard any residual complex parts
    all_ti_profiles(all_ti_profiles < 0) = 0;      % Temperature cannot be negative

    n_H0_profile = neutral_profile.n_H0(:); % Ensure neutral profile is a column vector

    % --- 2. Compute Charge Exchange Rate Coefficient ---
    % <sigma*v>_cx = T_i^0.318 * 1e-8 [cm^3/s], where T_i is in eV.
    % Formula from Cornelis et al. (1994), as cited in the thesis (Eq. 3.17).
    fprintf('... calculating charge exchange rate coefficient\n');
    all_cx_rates = (all_ti_profiles .^ 0.318) * 1e-8; % [cm^3/s]

    % Calculate statistics from the distribution of rate profiles (operating on rows)
    cx_rate_avg = mean(all_cx_rates, 1)';
    ci_percentiles_rate = prctile(all_cx_rates, [2.5, 97.5], 1);
    cx_rate_ci_lower = ci_percentiles_rate(1, :)';
    cx_rate_ci_upper = ci_percentiles_rate(2, :)';

    % --- 3. Compute Ion-Neutral Collision Frequency ---
    % nu_iH0 = <sigma*v>_cx * n_H0
    % A unit conversion from cm^3/s to m^3/s is required (factor of 1e-6).
    fprintf('... calculating ion-neutral collision frequency\n');
    all_nu_iH0 = (all_cx_rates * 1e-6) .* n_H0_profile'; % Broadcasting (matrix .* row vector)
    
    nu_iH0_avg = mean(all_nu_iH0, 1)';
    ci_nu = prctile(all_nu_iH0, [2.5, 97.5], 1);
    nu_iH0_ci_lower = ci_nu(1, :)';
    nu_iH0_ci_upper = ci_nu(2, :)';

    % --- 4. Package Results ---
    collision_profiles.r_fine = r_fine(:);
    
    collision_profiles.cx_rate.avg      = cx_rate_avg;
    collision_profiles.cx_rate.ci_lower = cx_rate_ci_lower;
    collision_profiles.cx_rate.ci_upper = cx_rate_ci_upper;
    
    collision_profiles.nu_iH0.avg       = nu_iH0_avg;
    collision_profiles.nu_iH0.ci_lower  = nu_iH0_ci_lower;
    collision_profiles.nu_iH0.ci_upper  = nu_iH0_ci_upper;

    fprintf('Collision profiles calculated successfully.\n\n');
end