% =========================================================================
% TCABR MOMENTUM ANALYSIS - MAIN SCRIPT
% =========================================================================
% This script serves as the main orchestrator for the analysis of toroidal
% momentum transport in the TCABR tokamak. It executes all calculation
% stages sequentially and then calls a separate script to generate plots.
%
% It includes a caching system to avoid re-running expensive calculations.
% =========================================================================

clc; clear; close all;
fprintf('Initiating TCABR momentum analysis script...\n\n');

addpath('src', 'plotting', 'utils');

%% --- SCRIPT CONTROL ---
% Set this flag to 'true' to force the script to run all calculations from
% scratch, even if a results file exists. Set to 'false' to load existing
% results and skip calculations.
force_recalculation = false;

% Define the path for the saved results file
results_filepath = fullfile('results', 'full_analysis_results.mat');


%% --- EXECUTION STAGES ---
% Check if results should be loaded or recalculated
if isfile(results_filepath) && ~force_recalculation
    % --- LOAD PRE-COMPUTED RESULTS ---
    fprintf('Found existing results file. Loading pre-computed data...\n');
    load(results_filepath);
    fprintf('Results loaded successfully.\n\n');
    
else
    % --- RUN ALL CALCULATION STAGES ---
    fprintf('Running all calculation stages from scratch...\n\n');
    
    % --- STAGE 1: Setup and Data Loading ---
    constants = setup_constants();
    exp_data = load_experimental_data(constants);
    r_fine = exp_data.r_profiles;
    
    % --- STAGE 2: Profile Analyses ---
    velocity_results = analyze_velocity_profile(exp_data, constants, r_fine);
    temperature_results = analyze_temperature_profile(exp_data, constants, r_fine);
    
    % --- STAGE 3: Physics Profile Calculations ---
    magnetic_field = compute_magnetic_field(r_fine, constants);
    derived_profiles = compute_derived_profiles(exp_data, temperature_results, magnetic_field, constants, r_fine);
    
    % --- STAGE 4: Theoretical Model Calculations ---
    theoretical_models = compute_theoretical_models(temperature_results, magnetic_field, derived_profiles, constants);
    
    % --- STAGE 5: Neutral Density and Collision Calculations ---
    neutral_profile = compute_neutral_density_profile(r_fine, constants);
    collision_profiles = compute_collision_profiles(temperature_results, neutral_profile, r_fine);
    
    % --- STAGE 6: Effective Diffusivity Calculation (Thesis Method) ---
    effective_diffusivity_thesis = compute_effective_diffusivity(velocity_results, collision_profiles, r_fine);
    
    % --- STAGE 7: Effective Diffusivity Calculation (from Scaling Laws) ---
    scaling_law_results = compute_chi_eff_from_scalings(velocity_results, derived_profiles, constants, r_fine);
    
    fprintf('\nAll calculations complete.\n');
    
    % --- SAVE RESULTS TO FILE ---
    fprintf('Saving all results to %s for future runs...\n', results_filepath);
    save(results_filepath, 'constants', 'exp_data', 'r_fine', ...
        'velocity_results', 'temperature_results', 'magnetic_field', ...
        'derived_profiles', 'theoretical_models', 'neutral_profile', ...
        'collision_profiles', 'effective_diffusivity_thesis', 'scaling_law_results');
    fprintf('Save complete.\n\n');
end

% --- STAGE 8: Plotting ---
% This stage always runs, using either the newly calculated or loaded results.
plot_diffusivity_comparison(effective_diffusivity_thesis, scaling_law_results, constants, r_fine);
% run_verification_plots(...); % You can uncomment this to see all plots

fprintf('Script finished.\n');