function run_scan_neutrals_fit(velocity_results, temperature_results, collision_profiles, ...
    theoretical_variants, constants, r_fine, results_dir)
%RUN_SCAN_NEUTRALS_FIT Configures and executes the neutral density parameter scan.
%   This function defines the parameter space for the grid search and selects
%   the theoretical model to be used as the reference for the fit. It then
%   calls the core scanning engine.

    fprintf('--- Configuring neutral density parameter scan ---\n');

    % --- 1. Define Parameter Ranges for the Scan ---
    amps = linspace(1.5e16, 3.0e16, 10);  % Amplitude of Gaussian peak [m^-3]
    widths = linspace(0.018, 0.028, 10); % Gaussian width (penetration depth) [m]
    r0s = linspace(0.15, 0.165, 5);      % Radial position of peak [m]
    
    % --- 2. Configure Scan Options ---
    opts = struct();
    opts.plot_results = true; % Generate and save plots from the scan
    opts.verbose = true;      % Print progress to the console
    
    % --- 3. Select the Theoretical Reference for Comparison ---
    % Use the 'Model.Variant' naming convention from the variants struct.
    % This is the "ground truth" we are fitting our experimental curve to.
    opts.reference_theory = 'Solomon.R0_mid'; 

    % --- 4. Execute the Scan ---
    scan_results = scan_neutrals_fit(velocity_results, temperature_results, collision_profiles, ...
        theoretical_variants, constants, r_fine, amps, widths, r0s, opts);

    % --- 5. Save Results ---
    scan_filename = fullfile(results_dir, 'study_neutral_scan_results.mat');
    fprintf('Saving neutral scan results to: %s\n', scan_filename);
    save(scan_filename, 'scan_results');
end