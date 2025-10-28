function run_scan_neutrals_fit(velocity_results, temperature_results, theoretical_variants, ...
    constants, r_fine, results_dir)
%RUN_SCAN_NEUTRALS_FIT Configures and executes the neutral density parameter scan.
%   This function defines the parameter space for the grid search and selects
%   the theoretical model to be used as the reference for the fit. It then
%   calls the core scanning engine.

    arguments
        velocity_results (1,1) struct, temperature_results (1,1) struct
        theoretical_variants (1,1) struct, constants (1,1) struct
        r_fine (:,1) double, results_dir (1,1) string
    end

    fprintf('--- Configuring neutral density parameter scan ---\n');

    % --- 1. Define Parameter Ranges for the Scan ---
    amps = linspace(1.5e16, 3.0e16, 10);  % Amplitude of Gaussian peak [m^-3]
    widths = linspace(0.020, 0.030, 10); % Gaussian width (penetration depth) [m]
    r0s = linspace(0.160, 0.165, 5);      % Radial position of peak [m]
    
    % --- 2. Configure Scan Options ---
    opts = struct();
    opts.plot_results = true; % Generate and save plots from the scan
    opts.verbose = true;      % Print progress to the console
    
    % --- 3. Select the Theoretical Reference for Comparison ---
    opts.reference_theory = 'Solomon.R0_mid'; 

    % --- 4. Execute the Scan ---
    % Note: 'collision_profiles' is removed as it's not used.
    scan_results = scan_neutrals_fit(velocity_results, temperature_results, ...
        theoretical_variants, constants, r_fine, amps, widths, r0s, opts);

    % --- 5. Save Results ---
    scan_filename = fullfile(results_dir, 'study_neutral_scan_results.mat');
    fprintf('Saving neutral scan results to: %s\n', scan_filename);
    save(scan_filename, 'scan_results');
end