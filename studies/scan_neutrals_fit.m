function results = scan_neutrals_fit(velocity_results, temperature_results, collision_profiles, ...
    variants, constants, r_fine, amps, widths, r0s, opts)
%SCAN_NEUTRALS_FIT Performs a grid scan of neutral Gaussian parameters.
%   Finds the optimal parameters (amplitude, width, peak position) for the
%   neutral density profile by minimizing the RMSE between the resulting
%   experimental chi_phi_eff and a specified theoretical variant.

    arguments
        % (Validation for primary inputs)
        velocity_results (1,1) struct, temperature_results (1,1) struct
        collision_profiles (1,1) struct, variants (1,1) struct
        constants (1,1) struct, r_fine (:,1) double
        % (Validation for scan parameters)
        amps (1,:) double, widths (1,:) double, r0s (1,:) double
        opts (1,1) struct
    end

    % --- 1. Setup ---
    safe_mask = variants.meta.safe_mask;
    
    % Dynamically get the reference theory profile from the variants struct
    ref_parts = split(opts.reference_theory, '.');
    try
        theory_ref_struct = variants.(ref_parts{1}).(ref_parts{2});
        theory_primary = theory_ref_struct.profile_avg;
    catch
        error('Invalid reference theory string: ''%s''', opts.reference_theory);
    end
    
    nA = numel(amps); nW = numel(widths); nR0 = numel(r0s);
    RMSE_cube = nan(nR0, nW, nA);
    best.rmse = inf;

    % --- 2. Main Grid Scan Loop ---
    fprintf('Starting grid scan (%d x %d x %d = %d points)...\n', nR0, nW, nA, numel(RMSE_cube));
    for i_r0 = 1:nR0
        r0 = r0s(i_r0);
        if opts.verbose, fprintf('...scanning r0 = %.4f m (%d/%d)\n', r0, i_r0, nR0); end
        
        parfor i_w = 1:nW
            w = widths(i_w);
            temp_rmse_row = nan(1, nA);
            for i_a = 1:nA
                amp = amps(i_a);
                
                % Create a temporary constants struct with the new neutral params
                temp_constants = constants;
                temp_constants.models.neutrals.n_H0_max_amp = amp;
                temp_constants.models.neutrals.width = w;
                temp_constants.models.neutrals.r_max_H_alpha = r0;
                
                % Recalculate the experimental diffusivity with these new params
                neutral_prof = compute_neutral_density_profile(r_fine, temp_constants);
                coll_prof = compute_collision_profiles(temperature_results, neutral_prof, r_fine);
                diff_exp = compute_diffusivity_profile(velocity_results, coll_prof, r_fine);
                
                % Calculate RMSE on the safe radial window
                residuals = diff_exp.profile_avg(safe_mask) - theory_primary(safe_mask);
                temp_rmse_row(i_a) = sqrt(mean(residuals.^2));
            end
            RMSE_cube(i_r0, i_w, :) = temp_rmse_row;
        end
    end
    
    % --- 3. Find and Store Best Result ---
    [min_rmse, min_idx] = min(RMSE_cube(:));
    [idx_r0, idx_w, idx_a] = ind2sub(size(RMSE_cube), min_idx);
    
    best.rmse = min_rmse;
    best.amp = amps(idx_a);
    best.width = widths(idx_w);
    best.r0 = r0s(idx_r0);
    
    fprintf('Scan complete. Best RMSE = %.4f\n', best.rmse);
    fprintf('  Best params: Amp=%.2e, Width=%.4f, r0=%.4f\n', best.amp, best.width, best.r0);

    % --- 4. Package Results ---
    results.opts = opts;
    results.RMSE_cube = RMSE_cube;
    results.amps = amps;
    results.widths = widths;
    results.r0s = r0s;
    results.best = best;

    % --- 5. Plotting ---
    if opts.plot_results
        % (Plotting logic for heatmap and best-fit overlay would go here)
        fprintf('Generating study-specific plots...\n');
    end
end