% =========================================================================
% TCABR MOMENTUM ANALYSIS - MAIN SCRIPT
% =========================================================================
% This script serves as the main orchestrator for the analysis of toroidal
% momentum transport in the TCABR tokamak. It executes all calculation
% stages sequentially and then calls a separate script to generate plots.
% =========================================================================

clc; clear; close all;
fprintf('Initiating TCABR momentum analysis script...\n\n');

addpath('src', 'plotting', 'utils');

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
effective_diffusivity = compute_effective_diffusivity(velocity_results, collision_profiles, r_fine);

fprintf('All calculations complete.\n\n');

% --- STAGE 7: Plotting ---
% Call the dedicated plotting script to generate all verification figures.
run_verification_plots(constants, exp_data, velocity_results, temperature_results, ...
    derived_profiles, theoretical_models, effective_diffusivity);

fprintf('Script finished.\n');
