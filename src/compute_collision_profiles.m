function [collision_profiles] = compute_collision_profiles(temperature_results, neutral_profile, r_fine)
%COMPUTE_COLLISION_PROFILES Calculates charge exchange rate and ion-neutral collision frequency.
%   This function robustly computes the charge exchange (cx) rate coefficient
%   and the ion-neutral collision frequency (nu_iH0). Crucially, it propagates
%   the uncertainty from the bootstrapped ion temperature profiles to generate
%   a 95% confidence interval for the output profiles.
%
%   Syntax:
%       collision_profiles = compute_collision_profiles(temperature_results, neutral_profile, r_fine)

    arguments
        temperature_results (1,1) struct {mustContainFields(temperature_results, {'bootstrap'})}
        neutral_profile (1,1) struct {mustContainFields(neutral_profile, {'n_H0'})}
        r_fine (:,1) {mustBeNumeric, mustBeReal, mustBeFinite}
    end

    fprintf('Calculating collision rate and frequency profiles...\n');

    % --- 1. Input Validation and Data Sanitisation ---
    mustContainFields(temperature_results.bootstrap, {'all_fitted_profiles'});
    all_ti_profiles = temperature_results.bootstrap.all_fitted_profiles;

    validateattributes(all_ti_profiles, {'numeric'}, {'2d', 'real', 'finite', 'nonnegative'}, mfilename, 'all_fitted_profiles');
    if size(all_ti_profiles, 2) ~= numel(r_fine)
        error('Dimension mismatch: The number of points in temperature profiles (%d) does not match the radial grid (%d).', ...
            size(all_ti_profiles, 2), numel(r_fine));
    end

    n_H0_profile = neutral_profile.n_H0(:);

    % --- 2. Compute Charge Exchange Rate Coefficient ---
    % Formula from Cornelis et al. (1994), as cited in the thesis (Eq. 3.17).
    % <sigma*v>_cx = T_i^0.318 * 1e-8 [cm^3/s], where T_i is in eV.
    fprintf('... calculating charge exchange rate coefficient\n');
    all_cx_rates = (all_ti_profiles .^ 0.318) * 1e-8; % [cm^3/s]

    % --- 3. Compute Ion-Neutral Collision Frequency ---
    % nu_iH0 = <sigma*v>_cx * n_H0, with unit conversion from cm^3 to m^3
    fprintf('... calculating ion-neutral collision frequency\n');
    all_nu_iH0 = (all_cx_rates * 1e-6) .* n_H0_profile'; % Broadcasting

    % --- 4. Calculate Statistics and Package Results ---
    collision_profiles.r_fine = r_fine;
    
    % CX Rate Statistics
    collision_profiles.cx_rate.avg = mean(all_cx_rates, 1)';
    ci_percentiles_rate = prctile(all_cx_rates, [2.5, 97.5], 1);
    collision_profiles.cx_rate.ci_lower = ci_percentiles_rate(1, :)';
    collision_profiles.cx_rate.ci_upper = ci_percentiles_rate(2, :)';
    
    % nu_iH0 Statistics
    collision_profiles.nu_iH0.avg = mean(all_nu_iH0, 1)';
    ci_percentiles_nu = prctile(all_nu_iH0, [2.5, 97.5], 1);
    collision_profiles.nu_iH0.ci_lower = ci_percentiles_nu(1, :)';
    collision_profiles.nu_iH0.ci_upper = ci_percentiles_nu(2, :)';

    fprintf('Collision profiles calculated successfully.\n\n');
end