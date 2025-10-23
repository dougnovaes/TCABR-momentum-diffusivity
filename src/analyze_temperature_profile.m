function [temp_results] = analyze_temperature_profile(exp_data, constants, r_fine)
%ANALYZE_TEMPERATURE_PROFILE Performs a detailed statistical analysis of the ion temperature profile.
%   This function fits a canonical profile model of the form:
%       T_i(r) = (T0 - Ta) * (1 - (r/a)^2)^sigma + Ta
%   to the experimental ion temperature data. It uses a robust Monte Carlo
%   bootstrap method for rigorous uncertainty quantification.
%
%   Key tasks:
%   1.  Performs a primary weighted non-linear least-squares fit to find the
%       optimal parameters (T0, Ta, sigma).
%   2.  Executes a bootstrap analysis by repeatedly fitting the model to data
%       perturbed within their error bars.
%   3.  Calculates the mean fitted profile, R-squared value, 95%
%       confidence band, and the uncertainty on the fitted parameters.
%
%   Syntax:
%       temp_results = analyze_temperature_profile(exp_data, constants, r_fine)

    arguments
        exp_data (1,1) struct
        constants (1,1) struct
        r_fine (:,1) {mustBeNumeric, mustBeReal, mustBeFinite}
    end

    fprintf('Starting ion temperature profile analysis...\n');
    t_start = tic;

    % --- 1. Input Validation and Preparation ---
    r_exp        = exp_data.r_Ti_exp(:);
    ti_recon_exp = exp_data.Ti_reconstructed_exp(:);
    ti_recon_err = exp_data.Ti_reconstructed_err_exp(:);
    
    validateattributes(r_exp, {'numeric'}, {'vector', 'real', 'finite', 'nonnegative'}, mfilename, 'exp_data.r_Ti_exp');
    validateattributes(ti_recon_exp, {'numeric'}, {'vector', 'real', 'finite', 'size', size(r_exp)}, mfilename, 'exp_data.Ti_reconstructed_exp');
    validateattributes(ti_recon_err, {'numeric'}, {'vector', 'real', 'finite', 'nonnegative', 'size', size(r_exp)}, mfilename, 'exp_data.Ti_reconstructed_err_exp');

    a            = constants.machine.a;
    n_iter       = constants.analysis.num_iterations;
    n_data       = numel(r_exp);

    % --- 2. Define Fit Model, Guesses, and Solver Options ---
    Ti_function = @(params, r_val) real_profile(params, r_val, a);

    initial_guess = [270; 30; 3.44]; % [T0; Ta; sigma]
    lower_bounds  = [-Inf; -Inf; 0];
    upper_bounds  = [Inf;  Inf;  Inf];

    weights_sqrt = 1 ./ ti_recon_err;
    weighted_residual = @(params) (Ti_function(params, r_exp) - ti_recon_exp) .* weights_sqrt;

    options = optimset('Display', 'off', 'TolFun', 1e-6, 'TolX', 1e-6);

    % --- 3. Perform Primary Fit on Unperturbed Data ---
    fprintf('Performing primary fit for optimal parameters...\n');
    try
        params_optimized = lsqnonlin(weighted_residual, initial_guess, lower_bounds, upper_bounds, options);
    catch ME
        % The following syntax for 'warning' is correct. The linter may be flagging a false positive.
        warning(ME.identifier, 'Primary fit failed: %s. Using initial guess as optimised parameters.', ME.message);
        params_optimized = initial_guess;
    end
    ti_fit_primary = Ti_function(params_optimized, r_fine);

    % Calculate R-squared for the primary fit
    ti_fit_at_data_points = Ti_function(params_optimized, r_exp);
    ss_res = sum(weights_sqrt.^2 .* (ti_recon_exp - ti_fit_at_data_points).^2);
    ss_tot = sum(weights_sqrt.^2 .* (ti_recon_exp - mean(ti_recon_exp)).^2);
    r_squared = 1 - ss_res / (ss_tot + eps);
    fprintf('Primary fit complete. R-squared: %.4f\n', r_squared);

    % --- 4. Monte Carlo Bootstrap for Uncertainty Analysis ---
    fprintf('Running %d bootstrap iterations for uncertainty analysis...\n', n_iter);
    ti_perturbed_sets = ti_recon_exp + ti_recon_err .* randn(n_data, n_iter);

    all_fitted_params   = nan(n_iter, numel(initial_guess));
    all_fitted_profiles = nan(n_iter, numel(r_fine));
    convergence_flags   = false(n_iter, 1);

    parfor i = 1:n_iter
        current_residual_func = @(params) ...
            (Ti_function(params, r_exp) - ti_perturbed_sets(:, i)) .* weights_sqrt;
        
        try
            params_perturbed = lsqnonlin(current_residual_func, initial_guess, lower_bounds, upper_bounds, options);
            if all(isfinite(params_perturbed))
                all_fitted_params(i, :) = params_perturbed;
                all_fitted_profiles(i, :) = Ti_function(params_perturbed, r_fine)';
                convergence_flags(i) = true;
            end
        catch
            % If lsqnonlin fails, pre-allocated NaNs are kept.
        end
    end

    % --- 5. Post-Process Bootstrap Results ---
    n_converged = sum(convergence_flags);
    % The following 'fprintf' syntax is correct. The linter may be flagging a false positive.
    fprintf('Bootstrap finished: %d of %d iterations succeeded (%.1f%%).\n', n_converged, n_iter, 100*n_converged/n_iter);
    if n_converged < 0.5 * n_iter
        warning('Less than 50%% of bootstrap iterations converged. Results may be unreliable.');
    elseif n_converged == 0
        error('All bootstrap iterations failed. Check model, initial guess, and data quality.');
    end

    % Filter valid results for both profiles and parameters
    valid_profiles = all_fitted_profiles(convergence_flags, :);
    valid_params   = all_fitted_params(convergence_flags, :);

    % Calculate profile statistics
    ti_profile_avg = mean(valid_profiles, 1)';
    ci_percentiles = prctile(valid_profiles, [2.5, 97.5], 1);
    ti_profile_ci_lower = ci_percentiles(1, :)';
    ti_profile_ci_upper = ci_percentiles(2, :)';
    
    % **CORRECTION**: Calculate parameter uncertainties from valid fits
    param_uncertainties = std(valid_params, 0, 1);

    % --- 6. Package Results ---
    temp_results.r_fine = r_fine;
    temp_results.profile_primary_fit = ti_fit_primary;
    temp_results.params_optimized.T0 = params_optimized(1);
    temp_results.params_optimized.Ta = params_optimized(2);
    temp_results.params_optimized.sigma = params_optimized(3);
    temp_results.r_squared = r_squared;
    temp_results.profile_avg = ti_profile_avg;
    temp_results.profile_ci_lower = ti_profile_ci_lower;
    temp_results.profile_ci_upper = ti_profile_ci_upper;
    
    % **CORRECTION**: Store the calculated parameter uncertainties
    temp_results.param_uncertainties.T0_err = param_uncertainties(1);
    temp_results.param_uncertainties.Ta_err = param_uncertainties(2);
    temp_results.param_uncertainties.sigma_err = param_uncertainties(3);
    
    temp_results.bootstrap.all_fitted_profiles = valid_profiles;
    temp_results.bootstrap.all_fitted_params = valid_params; % Storing for diagnostics
    temp_results.bootstrap.n_converged = n_converged;

    fprintf('Ion temperature analysis finished in %.2f s.\n\n', toc(t_start));
end

function Ti = real_profile(params, r, a)
    %REAL_PROFILE Safely evaluates the canonical Ti profile.
    T0 = params(1); Ta = params(2); sigma = params(3);
    base = 1 - (r(:) ./ a).^2;
    base(base < 0) = 0; % Prevents complex results for non-integer sigma
    Ti = (T0 - Ta) .* (base .^ sigma) + Ta;
    Ti = real(double(Ti(:)));
end 