function [velocity_results] = analyze_velocity_profile(exp_data, constants)
%ANALYZE_VELOCITY_PROFILE Performs a detailed statistical analysis of the toroidal velocity profile.
%   This function executes a robust analysis of the experimental toroidal
%   velocity data. Its main tasks are:
%   1.  Performs a Monte Carlo bootstrap analysis by repeatedly fitting a 5th-order 
%       weighted polynomial to data points perturbed within their error bars.
%   2.  Calculates the mean polynomial fit, the R-squared value, and the 95% 
%       confidence interval for the velocity profile.
%   3.  Fits the mean velocity profile to a Fourier-Bessel series to determine
%       the eigenvalues (lambda_j) of the system.
%   4.  Calculates the coefficient for the effective momentum diffusivity based on
%       the derived eigenvalues.
%
%   Syntax:
%       velocity_results = analyze_velocity_profile(exp_data, constants)
%
%   Inputs:
%       exp_data  - Structure containing the experimental data points (r_Vphi_exp, Vphi_exp, etc.).
%       constants - Structure containing machine parameters and analysis settings.
%
%   Output:
%       velocity_results - A structure containing all analysis results, including fitted
%                          profiles, confidence intervals, and derived physics coefficients.

fprintf('Starting toroidal velocity profile analysis...\n');
tic;

% --- 1. Initial Setup and Parameter Extraction ---
r_exp       = exp_data.r_Vphi_exp;
v_phi_exp   = exp_data.Vphi_exp * 1000;   % Convert km/s to m/s for calculations
v_phi_err   = exp_data.Vphi_err_exp * 1000; % Convert km/s to m/s
n_iter      = constants.analysis.num_iterations;
a           = constants.machine.a;
epsilon     = constants.machine.epsilon_aspect_ratio;

% Define a fine radial grid for high-resolution profile results
r_fine = linspace(eps, a, constants.analysis.num_radial_points)';

% --- 2. Monte Carlo Bootstrap for Polynomial Fit ---
% This loop fits a weighted polynomial for n_iter random perturbations
% of the experimental data to determine the uncertainty of the fit.

% Pre-allocate memory for performance
poly_coeffs_store = zeros(n_iter, 6); % 6 coeffs for a 5th-order polynomial

% Set up the fit type and options for the curve fitting toolbox
fit_type = fittype('poly5');
weights = 1 ./ (v_phi_err.^2);
fit_opts = fitoptions('Method', 'LinearLeastSquares', 'Weights', weights);

% Ensure a parallel pool is available for faster computation
if isempty(gcp('nocreate'))
    parpool;
end

fprintf('Running %d bootstrap iterations for polynomial fit...\n', n_iter);
parfor i = 1:n_iter
    % Perturb data points using a normal distribution defined by the error bars
    v_perturbed = v_phi_exp + v_phi_err .* randn(size(v_phi_exp));
    
    % Perform the weighted fit using the modern 'fit' function
    fitted_model = fit(r_exp, v_perturbed, fit_type, fit_opts);
    
    % Store the polynomial coefficients
    poly_coeffs_store(i, :) = coeffvalues(fitted_model);
end

% --- 3. Process Polynomial Fit Results ---
% Calculate the average fit and confidence intervals from the bootstrap results.

% Average polynomial coefficients and the resulting mean velocity profile
mean_poly_coeffs = mean(poly_coeffs_store, 1);
v_phi_poly_avg_fine = polyval(mean_poly_coeffs, r_fine);

% Evaluate all fitted polynomials on the fine grid to get profile statistics
v_phi_perturbed_fine = zeros(length(r_fine), n_iter);
for i = 1:n_iter
    v_phi_perturbed_fine(:, i) = polyval(poly_coeffs_store(i, :), r_fine);
end

% Calculate 95% confidence intervals using the percentile method (robust)
confidence_level = 0.95;
lower_percentile = (1 - confidence_level) / 2 * 100;
upper_percentile = (1 + confidence_level) / 2 * 100;

v_phi_poly_ci_lower = prctile(v_phi_perturbed_fine, lower_percentile, 2);
v_phi_poly_ci_upper = prctile(v_phi_perturbed_fine, upper_percentile, 2);

% Calculate R-squared value for the mean fit against the original data
predicted_values_exp = polyval(mean_poly_coeffs, r_exp);
ss_res = sum(weights .* (v_phi_exp - predicted_values_exp).^2);
ss_tot = sum(weights .* (v_phi_exp - mean(v_phi_exp)).^2);
r_squared = 1 - (ss_res / ss_tot);

fprintf('Polynomial analysis complete. R-squared: %.2f\n', r_squared);

% --- 4. Fourier-Bessel Series Analysis ---
% Fit the mean velocity profile to a series of Bessel functions to find eigenvalues.
fprintf('Performing Fourier-Bessel series fit...\n');

num_zeros = constants.analysis.num_zeros_bessel;
lambda_zeros = findBesselZeros(num_zeros, 0) / a; % Calculate lambda_j = z_j / a

% Formulate the linear system A*x = b for the Bessel fit
A_bessel = zeros(length(r_fine), num_zeros);
for k = 1:num_zeros
    A_bessel(:, k) = besselj(0, lambda_zeros(k) * r_fine);
end

% Solve for the Bessel coefficients using linear least squares (A \ b)
bessel_coeffs = A_bessel \ v_phi_poly_avg_fine;

% Reconstruct the velocity profile from the Bessel series
v_phi_bessel_fit_fine = A_bessel * bessel_coeffs;

% --- 5. Calculate Physics Parameters ---
% Use the derived eigenvalues to calculate the diffusivity coefficient.
sum_inv_lambda_sq = sum(1 ./ lambda_zeros.^2);

% This is the coefficient from Eq. 3.13, which multiplies nu_iH0.
% The nu_iH0 profile itself will be calculated in a separate function.
chi_phi_coeff = (sum_inv_lambda_sq) * ((3 / (4 * epsilon)) - 1);

fprintf('Bessel analysis complete. Sum of inverse squared lambdas: %.4f m^2\n', sum_inv_lambda_sq);

% --- 6. Package Results into Output Structure ---
velocity_results.r_fine = r_fine; % Fine radial grid [m]

% Polynomial fit results (in m/s)
velocity_results.poly_fit_avg = v_phi_poly_avg_fine;
velocity_results.poly_fit_ci_lower = v_phi_poly_ci_lower;
velocity_results.poly_fit_ci_upper = v_phi_poly_ci_upper;
velocity_results.r_squared = r_squared;
velocity_results.poly_coeffs = mean_poly_coeffs;

% Bessel fit results (in m/s)
velocity_results.bessel_fit = v_phi_bessel_fit_fine;
velocity_results.lambda_zeros = lambda_zeros;
velocity_results.sum_inv_lambda_sq = sum_inv_lambda_sq;

% Derived physics parameter coefficients
velocity_results.chi_phi_coeff = chi_phi_coeff; % [m^2]

analysis_time = toc;
fprintf('Toroidal velocity analysis finished in %.2f seconds.\n\n', analysis_time);

end