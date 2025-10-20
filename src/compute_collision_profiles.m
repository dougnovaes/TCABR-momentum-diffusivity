function [collision_profiles] = compute_collision_profiles(temperature_results, neutral_profile, r_fine)
%COMPUTE_COLLISION_PROFILES Calculates charge exchange rate and collision frequency.
%   This function computes two key profiles for momentum transport analysis:
%   1.  The charge exchange (cx) rate coefficient, <sigma*v>_cx, based on an
%       empirical power-law fit to the ion temperature.
%   2.  The ion-neutral collision frequency, nu_iH0, which is the product of
%       the rate coefficient and the neutral density profile.
%
%   The function robustly propagates uncertainties from the ion temperature
%   profile through to both calculated profiles.
%
%   Syntax:
%       collision_profiles = compute_collision_profiles(temperature_results, neutral_profile, r_fine)
%
%   Inputs:
%       temperature_results - Structure from analyze_temperature_profile,
%                             containing all bootstrap-fitted Ti profiles.
%       neutral_profile     - Structure with the calculated n_H0 profile.
%       r_fine              - The high-resolution radial grid [m].
%
%   Output:
%       collision_profiles - A structure containing the mean profiles and
%                            confidence bands for the rate coefficient and
%                            collision frequency.

fprintf('Calculating collision rate and frequency profiles...\n');

% --- 1. Extract necessary data ---
% The bootstrap analysis for Ti produced a matrix of (iterations x points).
% We will use this to propagate the uncertainty.
all_ti_profiles = temperature_results.bootstrap.all_fitted_profiles;
n_H0_profile    = neutral_profile.n_H0(:)'; % Ensure n_H0 is a row vector for broadcasting

% --- 2. Calculate Charge Exchange Rate Coefficient ---
% This calculation uses the empirical formula from Cornelis et al. (1994),
% as cited in the thesis (Eq. 3.17).
% <sigma*v>_cx = T_i^0.318 * 1e-8 [cm^3/s], where T_i is in eV.
fprintf('... calculating charge exchange rate coefficient\n');

all_cx_rates = (all_ti_profiles .^ 0.318) * 1e-8; % [cm^3/s]

% Calculate statistics from the distribution of rate profiles
cx_rate_avg = mean(all_cx_rates, 1)'; % Mean profile, transposed to column

% Calculate 95% confidence intervals using the percentile method
ci_percentiles_rate = prctile(all_cx_rates, [2.5, 97.5], 1);
cx_rate_ci_lower = ci_percentiles_rate(1, :)';
cx_rate_ci_upper = ci_percentiles_rate(2, :)';

% --- 3. Calculate Ion-Neutral Collision Frequency ---
% nu_iH0 = <sigma*v>_cx * n_H0
% A unit conversion from cm^3/s to m^3/s is required (factor of 1e-6).
fprintf('... calculating ion-neutral collision frequency\n');

% Use broadcasting (row vector .* matrix) to get all frequency profiles
all_nu_iH0 = (all_cx_rates * 1e-6) .* n_H0_profile; % [s^-1]

% Calculate statistics from the distribution of frequency profiles
nu_iH0_avg = mean(all_nu_iH0, 1)'; % Mean profile, transposed to column

ci_percentiles_nu = prctile(all_nu_iH0, [2.5, 97.5], 1);
nu_iH0_ci_lower = ci_percentiles_nu(1, :)';
nu_iH0_ci_upper = ci_percentiles_nu(2, :)';

% --- 4. Package Results into Output Structure ---
collision_profiles.r_fine = r_fine;
% Rate coefficient results
collision_profiles.cx_rate.avg = cx_rate_avg;
collision_profiles.cx_rate.ci_lower = cx_rate_ci_lower;
collision_profiles.cx_rate.ci_upper = cx_rate_ci_upper;
% Collision frequency results
collision_profiles.nu_iH0.avg = nu_iH0_avg;
collision_profiles.nu_iH0.ci_lower = nu_iH0_ci_lower;
collision_profiles.nu_iH0.ci_upper = nu_iH0_ci_upper;

fprintf('Collision profiles calculated successfully.\n\n');

end