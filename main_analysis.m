% =========================================================================
% TCABR MOMENTUM TRANSPORT ANALYSIS PIPELINE
% =========================================================================
% AUTHOR: [Seu Nome/Grupo]
% DATE: Jan 2025
%
% DESCRIPTION:
%   This script orchestrates the complete analysis of toroidal momentum 
%   transport in the TCABR tokamak. It processes experimental data to 
%   determine the effective momentum diffusivity and benchmarks it against 
%   multiple theoretical scaling laws.
%
% ARCHITECTURE:
%   The pipeline is divided into modular STAGES (1 to 7). Each stage:
%   1. Executes a specific set of physics calculations via functions in 'src/'.
%   2. Caches the results in 'results/*.mat' files to avoid redundant re-computation.
%   3. Can be enabled/disabled via the 'run_flags' structure.
%
%   After computation, a separate PLOTTING section generates publication-
%   quality figures if requested.
%
% KEY METHODOLOGY (HYBRID APPROACH):
%   The theoretical benchmarking (Stage 7) uses a hybrid global/local method:
%   - Magnitude: Determined by global volume-averaged parameters (Maslov).
%   - Shape: Modulated by local experimental velocity gradients.
%
% REQUIREMENTS:
%   MATLAB R2021b+ (Curve Fitting Toolbox, Parallel Computing Toolbox).
% =========================================================================

clc; clear; close all;
fprintf('==============================================================\n');
fprintf('   TCABR MODULAR MOMENTUM ANALYSIS PIPELINE\n');
fprintf('==============================================================\n\n');

% -------------------------------------------------------------------------
% 1. ENVIRONMENT SETUP
% -------------------------------------------------------------------------
[base_dir, ~, ~] = fileparts(mfilename('fullpath'));

% Define module paths
addpath(fullfile(base_dir, 'src'));      % Core physics calculations
addpath(fullfile(base_dir, 'plotting')); % Visualization scripts
addpath(fullfile(base_dir, 'utils'));    % Helper utilities (cache, validation)

% Define output directories
results_dir = fullfile(base_dir, 'results');
plots_dir   = fullfile(results_dir, 'plots');

% Create directories if they don't exist
if ~exist(results_dir, 'dir'), mkdir(results_dir); end
if ~exist(plots_dir, 'dir'),   mkdir(plots_dir);   end

% -------------------------------------------------------------------------
% 2. CONTROL PANEL (USER CONFIGURATION)
% -------------------------------------------------------------------------
% Set 'true' to force re-calculation of a stage.
% Set 'false' to load results from the cache (.mat) if available.

run_flags = struct( ...
    'setup_and_data',        true, ...  % Stage 1: Load constants & raw data
    'profile_analyses',      true, ...  % Stage 2: Fit Vphi/Ti (Bootstrap)
    'physics_profiles',      true, ...  % Stage 3: Calc Global Collisionality & Local q, s
    'helander_model',        true, ...  % Stage 4: Calc Helander Velocity (Validation)
    'neutral_and_collision', true, ...  % Stage 5: Calc Neutral Density & Frequencies
    'diffusivity_thesis',    true, ...  % Stage 6: Calc Exp. Effective Diffusivity
    'diffusivity_theory',    true  ...  % Stage 7: Calc Theoretical Benchmarks (Hybrid)
);

% Visualization Controls
plot_flags = struct( ...
    'justification_and_fits',  true, ... % Exp. Fits & Helander Validation
    'thesis_method_results',   true, ... % Neutrals, Freqs & Exp. Diffusivity
    'supporting_profiles',     true, ... % q, s, Thermal Velocities, Collisions
    'final_comparison',        true, ... % Main Results: Exp vs Theory (2 Figures)
    'collisionality_comparison', true ... % Diagnostic: Wesson vs Solomon
);

% Execution Behavior
stop_on_error = true; % Halt pipeline immediately if a stage fails

fprintf('>>> Configuration loaded.\n\n');

% =========================================================================
% 3. COMPUTATION PIPELINE
% =========================================================================

%% STAGE 1: SETUP & DATA LOADING
%  Loads physical constants and raw experimental profiles from text files.
try
    [constants, exp_data, r_fine] = load_or_compute( ...
        fullfile(results_dir, 'stage1_setup_data.mat'), ...
        run_flags.setup_and_data, ...
        @stage1_setup_and_data, ...
        {'constants', 'exp_data', 'r_fine'});
    
    fprintf('[✓] Stage 1 completed: Setup & Data Loading.\n');
catch ME
    handle_stage_error(ME, 1, stop_on_error);
end

%% STAGE 2: PROFILE ANALYSIS
%  Fits smooth polynomials to V_phi and T_i data using Bootstrap for uncertainty.
try
    [velocity_results, temperature_results] = load_or_compute( ...
        fullfile(results_dir, 'stage2_profiles.mat'), ...
        run_flags.profile_analyses, ...
        @() stage2_profiles(exp_data, constants, r_fine), ...
        {'velocity_results', 'temperature_results'});
    
    fprintf('[✓] Stage 2 completed: Profile Analysis (Fits).\n');
catch ME
    handle_stage_error(ME, 2, stop_on_error);
end

%% STAGE 3: PHYSICS PROFILES
%  Calculates q, shear, thermal velocities, and CRUCIALLY the Global Maslov
%  Collisionality used for scaling laws.
try
    [magnetic_field, derived_profiles] = load_or_compute( ...
        fullfile(results_dir, 'stage3_physics.mat'), ...
        run_flags.physics_profiles, ...
        @() stage3_physics(exp_data, temperature_results, constants, r_fine), ...
        {'magnetic_field', 'derived_profiles'});
    
    fprintf('[✓] Stage 3 completed: Physics Profiles (Local & Global).\n');
catch ME
    handle_stage_error(ME, 3, stop_on_error);
end

%% STAGE 4: HELANDER MODEL VALIDATION
%  Calculates the neoclassical velocity prediction to validate the physics basis.
try
    helander_results = load_or_compute( ...
        fullfile(results_dir, 'stage4_helander.mat'), ...
        run_flags.helander_model, ...
        @() compute_helander_model(temperature_results, magnetic_field, derived_profiles, constants), ...
        {'helander_results'});
    
    fprintf('[✓] Stage 4 completed: Helander Model Validation.\n');
catch ME
    handle_stage_error(ME, 4, stop_on_error);
end

%% STAGE 5: NEUTRALS & COLLISIONS
%  Calculates the neutral density profile and ion-neutral collision frequencies.
try
    [neutral_profile, collision_profiles] = load_or_compute( ...
        fullfile(results_dir, 'stage5_collisions.mat'), ...
        run_flags.neutral_and_collision, ...
        @() stage5_neutral_collision(r_fine, temperature_results, constants), ...
        {'neutral_profile', 'collision_profiles'});
    
    fprintf('[✓] Stage 5 completed: Neutrals & Collision Frequencies.\n');
catch ME
    handle_stage_error(ME, 5, stop_on_error);
end

%% STAGE 6: EXPERIMENTAL EFFECTIVE DIFFUSIVITY
%  Calculates chi_eff using the Thesis Method (Bessel Eigenvalues).
try
    effective_diffusivity_thesis = load_or_compute( ...
        fullfile(results_dir, 'stage6_diffusivity_thesis.mat'), ...
        run_flags.diffusivity_thesis, ...
        @() compute_diffusivity_profile(velocity_results, collision_profiles, r_fine), ...
        {'effective_diffusivity_thesis'});
    
    fprintf('[✓] Stage 6 completed: Experimental Diffusivity (Thesis).\n');
catch ME
    handle_stage_error(ME, 6, stop_on_error);
end

%% STAGE 7: THEORETICAL BENCHMARKING
%  Calculates theoretical transport profiles (Solomon, Hahm, Gürcan, Peeters)
%  using the Hybrid Global/Local approach.
try
    theoretical_variants = load_or_compute( ...
        fullfile(results_dir, 'stage7_theoretical_transport.mat'), ...
        run_flags.diffusivity_theory, ...
        @() compute_theoretical_transport(velocity_results, derived_profiles, constants, r_fine), ...
        {'theoretical_variants'});
    
    fprintf('[✓] Stage 7 completed: Theoretical Benchmarks.\n');
catch ME
    handle_stage_error(ME, 7, stop_on_error);
end

fprintf('\n>>> All computational stages executed successfully.\n\n');

% =========================================================================
% 4. VISUALIZATION SECTION
% =========================================================================
% Check if any plots are requested
p_vals = struct2cell(plot_flags);
if any([p_vals{:}])
    fprintf('--- Generating Selected Plots ---\n');
    
    try
        % 1. Justification & Fits (Exp vs Helander)
        if plot_flags.justification_and_fits
            plot_stage1_justification_and_fits(constants, exp_data, velocity_results, ...
                temperature_results, helander_results, r_fine, plots_dir);
        end

        % 2. Thesis Method Results (Neutrals & Exp. Chi)
        if plot_flags.thesis_method_results
            plot_stage2_thesis_method(constants, neutral_profile, collision_profiles, ...
                effective_diffusivity_thesis, r_fine, plots_dir);
        end

        % 3. Supporting Physics (q, s, nu*, etc.)
        if plot_flags.supporting_profiles
            plot_stage3_supporting_profiles(constants, derived_profiles, r_fine, plots_dir);
        end

        % 4. Final Comparison (Exp vs Theory - The Core Result)
        if plot_flags.final_comparison
            plot_final_comparison_variants(effective_diffusivity_thesis, theoretical_variants, ...
                constants, plots_dir);
        end
        
        % 5. Collisionality Sensitivity (Diagnostic)
        if plot_flags.collisionality_comparison
             plot_collisionality_comparison(derived_profiles, constants, velocity_results, ...
                effective_diffusivity_thesis, plots_dir);
        end

        fprintf('--- Plotting complete. Figures saved to: %s ---\n', plots_dir);
        
    catch ME
        warning('TCABR:PlottingError', 'An error occurred during plotting: %s', ME.message);
        fprintf(2, 'Error location: %s (line %d)\n', ME.stack(1).file, ME.stack(1).line);
    end
end

fprintf('\n>>> Pipeline Finished Successfully.\n');


% =========================================================================
% LOCAL WRAPPER FUNCTIONS (Data Passing Interfaces)
% =========================================================================

function [constants, exp_data, r_fine] = stage1_setup_and_data()
    constants = setup_constants();
    exp_data = load_experimental_data(constants);
    r_fine = exp_data.r_profiles;
end

function [velocity_results, temperature_results] = stage2_profiles(exp_data, constants, r_fine)
    velocity_results = analyze_velocity_profile(exp_data, constants, r_fine);
    temperature_results = analyze_temperature_profile(exp_data, constants, r_fine);
end

function [magnetic_field, derived_profiles] = stage3_physics(exp_data, temperature_results, constants, r_fine)
    magnetic_field = compute_magnetic_field(r_fine, constants);
    derived_profiles = compute_derived_profiles(exp_data, temperature_results, magnetic_field, constants, r_fine);
end

function [neutral_profile, collision_profiles] = stage5_neutral_collision(r_fine, temperature_results, constants)
    neutral_profile = compute_neutral_density_profile(r_fine, constants);
    collision_profiles = compute_collision_profiles(temperature_results, neutral_profile, r_fine);
end