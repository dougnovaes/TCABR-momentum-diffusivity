function [velocity_results] = analyze_velocity_profile(exp_data, constants, r_fine)
%ANALYZE_VELOCITY_PROFILE Performs a detailed statistical analysis of the toroidal velocity profile.
%   This function executes a two-stage analysis of the experimental V_phi data:
%   1.  **Polynomial Fit Bootstrap**: A 5th-order weighted polynomial is repeatedly
%       fit to the data, perturbed within error bars, to robustly determine the
%       mean velocity profile and its 95% confidence interval.
%   2.  **Fourier-Bessel Decomposition**: The mean velocity profile is fit to a
%       Fourier-Bessel series to determine the system's eigenvalues (lambda_j).
%       These are used to calculate the geometric coefficient of the effective
%       momentum diffusivity in the thesis model (Eq. 3.13).
%
%   Syntax:
%       velocity_results = analyze_velocity_profile(exp_data, constants, r_fine)
%
%   Inputs:
%       exp_data  - Structure with experimental Vphi data points.
%       constants - Structure with machine and analysis settings.
%       r_fine    - The high-resolution radial grid for profile evaluation [m].
%
%   Output:
%       velocity_results - A structure containing all analysis results.

    arguments
        exp_data (1,1) struct
        constants (1,1) struct
        r_fine (:,1) {mustBeNumeric, mustBeReal, mustBeFinite}
    end

    fprintf('Starting toroidal velocity profile analysis...\n');
    tic;

    % --- 1. Input Validation and Preparation ---
    r_exp       = exp_data.r_Vphi_exp(:);
    v_phi_exp   = exp_data.Vphi_exp(:) * 1000;   % Convert km/s to m/s
    v_phi_err   = exp_data.Vphi_err_exp(:) * 1000; % Convert km/s to m/s

    validateattributes(r_exp, {'numeric'}, {'vector', 'real', 'finite', 'nonnegative'}, mfilename, 'exp_data.r_Vphi_exp');
    validateattributes(v_phi_exp, {'numeric'}, {'vector', 'real', 'finite', 'size', size(r_exp)}, mfilename, 'exp_data.Vphi_exp');
    validateattributes(v_phi_err, {'numeric'}, {'vector', 'real', 'finite', 'nonnegative', 'size', size(r_exp)}, mfilename, 'exp_data.Vphi_err_exp');

    n_iter      = constants.analysis.num_iterations;
    a           = constants.machine.a;
    epsilon     = constants.machine.epsilon_aspect_ratio;

    % --- 2. Monte Carlo Bootstrap for Polynomial Fit ---
    fprintf('Running %d bootstrap iterations for polynomial fit...\n', n_iter);
    
    poly_coeffs_store = zeros(n_iter, 6); % 6 coeffs for a 5th-order polynomial
    fit_type = fittype('poly5');
    weights = 1 ./ (v_phi_err.^2);
    fit_opts = fitoptions('Method', 'LinearLeastSquares', 'Weights', weights);

    parfor i = 1:n_iter
        v_perturbed = v_phi_exp + v_phi_err .* randn(size(v_phi_exp));
        fitted_model = fit(r_exp, v_perturbed, fit_type, fit_opts);
        poly_coeffs_store(i, :) = coeffvalues(fitted_model);
    end

    % --- 3. Process Polynomial Fit Results ---
    mean_poly_coeffs = mean(poly_coeffs_store, 1);
    v_phi_poly_avg_fine = polyval(mean_poly_coeffs, r_fine);

    v_phi_all_fits_fine = zeros(numel(r_fine), n_iter);
    for i = 1:n_iter
        v_phi_all_fits_fine(:, i) = polyval(poly_coeffs_store(i, :), r_fine);
    end

    ci_percentiles = prctile(v_phi_all_fits_fine, [2.5, 97.5], 2);
    v_phi_poly_ci_lower = ci_percentiles(:, 1);
    v_phi_poly_ci_upper = ci_percentiles(:, 2);

    predicted_values_exp = polyval(mean_poly_coeffs, r_exp);
    ss_res = sum(weights .* (v_phi_exp - predicted_values_exp).^2);
    ss_tot = sum(weights .* (v_phi_exp - mean(v_phi_exp)).^2);
    r_squared = 1 - (ss_res / ss_tot);
    fprintf('Polynomial analysis complete. R-squared: %.4f\n', r_squared);

    % --- 4. Fourier-Bessel Series Analysis ---
    fprintf('Performing Fourier-Bessel series fit...\n');
    num_zeros = constants.analysis.num_zeros_bessel;
    lambda_zeros = findBesselZeros(num_zeros, 0) / a; % lambda_j = z_j / a

    A_bessel = zeros(length(r_fine), num_zeros);
    for k = 1:num_zeros
        A_bessel(:, k) = besselj(0, lambda_zeros(k) * r_fine);
    end

    bessel_coeffs = A_bessel \ v_phi_poly_avg_fine;
    v_phi_bessel_fit_fine = A_bessel * bessel_coeffs;

    % --- 5. Calculate Physics Parameters from Eigenvalues ---
    % This coefficient from Eq. 3.13 multiplies nu_iH0 in the thesis model.
    sum_inv_lambda_sq = sum(1 ./ lambda_zeros.^2);
    chi_phi_coeff = (sum_inv_lambda_sq) * ((3 / (4 * epsilon)) - 1);
    fprintf('Bessel analysis complete. chi_phi coefficient: %.4f m^2\n', chi_phi_coeff);

    % --- 6. Package Results ---
    velocity_results.r_fine = r_fine;
    velocity_results.poly_fit_avg = v_phi_poly_avg_fine;
    velocity_results.poly_fit_ci_lower = v_phi_poly_ci_lower;
    velocity_results.poly_fit_ci_upper = v_phi_poly_ci_upper;
    velocity_results.r_squared = r_squared;
    velocity_results.poly_coeffs = mean_poly_coeffs;
    velocity_results.bessel_fit = v_phi_bessel_fit_fine;
    velocity_results.lambda_zeros = lambda_zeros;
    velocity_results.chi_phi_coeff = chi_phi_coeff; % [m^2]

    fprintf('Toroidal velocity analysis finished in %.2f seconds.\n\n', toc);
end