function [temp_results] = analyze_temperature_profile(exp_data, constants)
%ANALYZE_TEMPERATURE_PROFILE Performs a detailed statistical analysis of the ion temperature profile.
%   This function fits a canonical profile model to the experimental ion
%   temperature data. Its main tasks are:
%   1.  Performs a primary weighted non-linear least-squares fit to the 
%       reconstructed experimental data to find the optimal parameters (T0, Ta, sigma).
%   2.  Executes a Monte Carlo bootstrap analysis by repeatedly fitting the model
%       to data points perturbed within their error bars to determine uncertainties.
%   3.  Calculates the mean fitted profile, R-squared value, and the 95% 
%       confidence band.
%
%   Syntax:
%       temp_results = analyze_temperature_profile(exp_data, constants)
%
%   Inputs:
%       exp_data  - Structure containing the experimental Ti data points.
%       constants - Structure containing machine parameters and analysis settings.
%
%   Output:
%       temp_results - A structure containing all analysis results for the ion temperature.

fprintf('Starting ion temperature profile analysis...\n');
tic;

% --- 1. Initial Setup and Parameter Extraction ---
r_exp        = exp_data.r_Ti_exp;
ti_recon_exp = exp_data.Ti_reconstructed_exp;
ti_recon_err = exp_data.Ti_reconstructed_err_exp;
a            = constants.machine.a;
n_iter       = constants.analysis.num_iterations;

% Define a fine radial grid for high-resolution profile results
r_fine = linspace(eps, a, constants.analysis.num_radial_points)';

% --- 2. Define Fit Model and Options ---
% The canonical profile function, as defined in the original script.
Ti_function = @(params, r_val) (params(1) - params(2)) .* (1 - (r_val / a).^2).^params(3) + params(2);

% Initial guess for the parameters [T0, Ta, sigma]
initial_guess = [270, 30, 3.44];

% Define the weighted objective function using 1/err^2 for reconstructed data
weights_sq = 1 ./ ti_recon_err.^2;
weighted_residual = @(params) (Ti_function(params, r_exp) - ti_recon_exp) .* sqrt(weights_sq);

% Set options for the non-linear solver
options = optimset('Display', 'off', 'TolFun', 1e-6, 'TolX', 1e-6);

% --- 3. Perform Primary Fit on Unperturbed Data ---
fprintf('Performing primary fit for optimal parameters...\n');
params_optimized = lsqnonlin(weighted_residual, initial_guess, [], [], options);
ti_fit_primary = Ti_function(params_optimized, r_fine);

% Calculate R-squared for the primary fit
ti_fit_at_data_points = Ti_function(params_optimized, r_exp);
ss_res = sum((ti_recon_exp - ti_fit_at_data_points).^2);
ss_tot = sum((ti_recon_exp - mean(ti_recon_exp)).^2);
r_squared = 1 - (ss_res / ss_tot);

fprintf('Primary fit complete. R-squared: %.3f\n', r_squared);

% --- 4. Monte Carlo Bootstrap for Uncertainty Analysis ---
% Generate all perturbed Ti data points in a vectorized way for efficiency
ti_perturbed_sets = ti_recon_exp + ti_recon_err .* randn(length(ti_recon_exp), n_iter);

% Pre-allocate memory
all_fitted_params = zeros(n_iter, length(initial_guess));
all_fitted_profiles = zeros(length(r_fine), n_iter);

% Ensure a parallel pool is available
if isempty(gcp('nocreate'))
    parpool;
end

fprintf('Running %d bootstrap iterations for uncertainty analysis...\n', n_iter);
parfor i = 1:n_iter
    % Define the weighted residual for the current perturbed data set
    current_residual_func = @(params) ...
        (Ti_function(params, r_exp) - ti_perturbed_sets(:, i)) .* sqrt(weights_sq);
    
    % Fit the perturbed data
    params_perturbed = lsqnonlin(current_residual_func, initial_guess, [], [], options);
    
    % Store the parameters and the full profile for this iteration
    all_fitted_params(i, :) = params_perturbed;
    all_fitted_profiles(:, i) = Ti_function(params_perturbed, r_fine);
end

% --- 5. Process Bootstrap Results ---
% Calculate mean profile, confidence intervals, and parameter uncertainties.

ti_profile_avg = mean(all_fitted_profiles, 2);
ti_profile_std = std(all_fitted_profiles, 0, 2);

% Calculate 95% confidence intervals using the percentile method
confidence_level = 0.95;
lower_percentile = (1 - confidence_level) / 2 * 100;
upper_percentile = (1 + confidence_level) / 2 * 100;

ti_profile_ci_lower = prctile(all_fitted_profiles, lower_percentile, 2);
ti_profile_ci_upper = prctile(all_fitted_profiles, upper_percentile, 2);

% Calculate the standard deviation of the parameters for uncertainty
param_uncertainties = std(all_fitted_params, 0, 1);

% --- 6. Package Results into Output Structure ---
temp_results.r_fine = r_fine; % Fine radial grid [m]

% Primary fit results
temp_results.profile_primary_fit = ti_fit_primary; % Profile from main fit [eV]
temp_results.params_optimized.T0 = params_optimized(1); % [eV]
temp_results.params_optimized.Ta = params_optimized(2); % [eV]
temp_results.params_optimized.sigma = params_optimized(3);
temp_results.r_squared = r_squared;

% Bootstrap analysis results
temp_results.profile_avg = ti_profile_avg; % Mean profile from all iterations [eV]
temp_results.profile_ci_lower = ti_profile_ci_lower; % [eV]
temp_results.profile_ci_upper = ti_profile_ci_upper; % [eV]
temp_results.param_uncertainties.T0_err = param_uncertainties(1); % [eV]
temp_results.param_uncertainties.Ta_err = param_uncertainties(2); % [eV]
temp_results.param_uncertainties.sigma_err = param_uncertainties(3);

analysis_time = toc;
fprintf('Ion temperature analysis finished in %.2f seconds.\n\n', analysis_time);

end