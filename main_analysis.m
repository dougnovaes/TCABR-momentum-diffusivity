% =========================================================================
% TCABR MOMENTUM TRANSPORT ANALYSIS (MODULAR & CACHE-ENABLED)
% =========================================================================
% PURPOSE:
%   This main script orchestrates the entire data analysis pipeline for
%   toroidal momentum transport in the TCABR tokamak. It is designed for
%   modular execution, per-stage caching, and controlled plotting.
%
% DESIGN:
%   - Modular "stage" system: Results are saved to 'results/' as .mat files.
%   - Focus: Calculation of effective diffusivity and comparison with 
%     theoretical benchmarks (Solomon, Hahm, Gürcan, Peeters).
%   - No inverse fitting or parameter scanning is performed in this version.
%
% LAYOUT:
%   Root/
%     ├── main_analysis.m    # This script
%     ├── data/              # Input data files (experimental_profiles.txt)
%     ├── src/               # Core computational routines
%     ├── utils/             # Utility functions (cache, math)
%     ├── plotting/          # Publication-quality visualisation
%     └── results/           # Output .mat files and figures
%
% REQUIREMENTS:
%   MATLAB R2025b or newer (Toolboxes: Curve Fitting, Parallel Computing)
% =========================================================================

clc; clear; close all;
fprintf('>>> Initiating TCABR modular momentum analysis pipeline...\n\n');

% -------------------------------------------------------------------------
% 1. PATH SETUP
% -------------------------------------------------------------------------
[base_dir, ~, ~] = fileparts(mfilename('fullpath'));
addpath(fullfile(base_dir, 'src'), ...
        fullfile(base_dir, 'plotting'), ...
        fullfile(base_dir, 'utils'));
        % Note: 'studies' folder removed from path as it is no longer used.

results_dir = fullfile(base_dir, 'results');
if ~exist(results_dir, 'dir'), mkdir(results_dir); end

plots_dir = fullfile(results_dir, 'plots');
if ~exist(plots_dir, 'dir'), mkdir(plots_dir); end

% -------------------------------------------------------------------------
% 2. CONTROL PANEL
% -------------------------------------------------------------------------
% CONFIGURATION:
% - Set flags to 'true' to force re-computation.
% - Set to 'false' to load from cache (if available).

run_flags = struct( ...
    'setup_and_data',        true, ... % Stage 1: Load constants & raw data
    'profile_analyses',      true, ... % Stage 2: Fit Vphi/Ti (Bootstrap)
    'physics_profiles',      true, ... % Stage 3: Calc Collisionality (Solomon Global), q, s
    'theoretical_models',    true, ... % Stage 4: Calc Theoretical Profiles (Helander/Solomon)
    'neutral_and_collision', true, ... % Stage 5: Calc Neutral Density & Freqs
    'diffusivity_thesis',    true, ... % Stage 6: Calc Exp. Chi_eff (Thesis Method)
    'diffusivity_scaling',   true, ... % Stage 7a: Calc Theory Chi (Base Scaling)
    'build_variants',        true   ... % Stage 7b: Calc Theory Variants (Global Magnitude/Local Shape)
);

% PLOTTING CONFIGURATION:
plot_flags = struct( ...
    'justification_and_fits',  true, ... % Stage 1 Plots
    'thesis_method_results',   true, ... % Stage 2 Plots
    'supporting_profiles',     true, ... % Stage 3 Plots
    'final_comparison',        true,  ... % Stage 4: Final Comparison (Fixed Scale 0-25)
    'collisionality_comparison', true ... % Analysis: Wesson vs Solomon
);

stop_on_error = true;

fprintf('>>> Control panel configured.\n\n');

% =========================================================================
% 3. PIPELINE EXECUTION (Calculations)
% =========================================================================

%% Stage 1: Setup & Data
try
    [constants, exp_data, r_fine] = load_or_compute( ...
        fullfile(results_dir, 'stage1_setup_data.mat'), ...
        run_flags.setup_and_data, ...
        @stage1_setup_and_data, ...
        {'constants','exp_data','r_fine'});
    fprintf('[✓] Stage 1 completed: Setup & Data.\n');
catch ME
    handle_stage_error(ME, 1, stop_on_error);
end

%% Stage 2: Profile Analyses
try
    [velocity_results, temperature_results] = load_or_compute( ...
        fullfile(results_dir, 'stage2_profiles.mat'), ...
        run_flags.profile_analyses, ...
        @() stage2_profiles(exp_data, constants, r_fine), ...
        {'velocity_results','temperature_results'});
    fprintf('[✓] Stage 2 completed: Profile Analyses.\n');
catch ME
    handle_stage_error(ME, 2, stop_on_error);
end

%% Stage 3: Physics Profiles
try
    [magnetic_field, derived_profiles] = load_or_compute( ...
        fullfile(results_dir, 'stage3_physics.mat'), ...
        run_flags.physics_profiles, ...
        @() stage3_physics(exp_data, temperature_results, constants, r_fine), ...
        {'magnetic_field','derived_profiles'});
    fprintf('[✓] Stage 3 completed: Physics Profiles.\n');
catch ME
    handle_stage_error(ME, 3, stop_on_error);
end

%% Stage 4: Theoretical Models
try
    theoretical_models = load_or_compute( ...
        fullfile(results_dir, 'stage4_models.mat'), ...
        run_flags.theoretical_models, ...
        @() stage4_models(temperature_results, magnetic_field, derived_profiles, constants), ...
        {'theoretical_models'});
    fprintf('[✓] Stage 4 completed: Theoretical Models.\n');
catch ME
    handle_stage_error(ME, 4, stop_on_error);
end

%% Stage 5: Neutrals & Collisions
try
    [neutral_profile, collision_profiles] = load_or_compute( ...
        fullfile(results_dir, 'stage5_collisions.mat'), ...
        run_flags.neutral_and_collision, ...
        @() stage5_neutral_collision(r_fine, temperature_results, constants), ...
        {'neutral_profile','collision_profiles'});
    fprintf('[✓] Stage 5 completed: Neutrals & Collisions.\n');
catch ME
    handle_stage_error(ME, 5, stop_on_error);
end

%% Stage 6: Diffusivity (Thesis Method)
try
    effective_diffusivity_thesis = load_or_compute( ...
        fullfile(results_dir, 'stage6_diffusivity_thesis.mat'), ...
        run_flags.diffusivity_thesis, ...
        @() stage6_diffusivity_thesis(velocity_results, collision_profiles, r_fine), ...
        {'effective_diffusivity_thesis'});
    fprintf('[✓] Stage 6 completed: Exp. Effective Diffusivity.\n');
catch ME
    handle_stage_error(ME, 6, stop_on_error);
end

%% Stage 7a: Diffusivity (Scaling Laws Base)
try
    scaling_law_results = load_or_compute( ...
        fullfile(results_dir, 'stage7a_diffusivity_scaling_base.mat'), ...
        run_flags.diffusivity_scaling, ...
        @() stage7a_diffusivity_scaling(velocity_results, derived_profiles, constants, r_fine), ...
        {'scaling_law_results'});
    fprintf('[✓] Stage 7a completed: Theor. Diffusivity (Base).\n');
catch ME
    handle_stage_error(ME, 7, stop_on_error);
end

%% Stage 7b: Build Theoretical Variants
% This stage calls the external function 'src/build_theoretical_variants.m'
try
    theoretical_variants = load_or_compute( ...
        fullfile(results_dir, 'stage7b_theoretical_variants.mat'), ...
        run_flags.build_variants, ...
        @() build_theoretical_variants(velocity_results, derived_profiles, constants, r_fine), ...
        {'theoretical_variants'});
    fprintf('[✓] Stage 7b completed: Theor. Variants Built.\n');
catch ME
    handle_stage_error(ME, 7.5, stop_on_error);
end

fprintf('\n>>> All calculation stages executed (or loaded) successfully.\n\n');

% =========================================================================
% 4. PLOTTING SECTION
% =========================================================================
p_cells = struct2cell(plot_flags);
p_vals = [p_cells{:}];

if any(p_vals)
    fprintf('--- Generating selected plots ---\n');
    try
        if plot_flags.justification_and_fits
            plot_stage1_justification_and_fits(constants, exp_data, velocity_results, ...
                temperature_results, theoretical_models, r_fine, plots_dir);
        end

        if plot_flags.thesis_method_results
            plot_stage2_thesis_method(constants, neutral_profile, collision_profiles, ...
                effective_diffusivity_thesis, r_fine, plots_dir);
        end

        if plot_flags.supporting_profiles
            plot_stage3_supporting_profiles(constants, derived_profiles, r_fine, plots_dir);
        end

        if plot_flags.final_comparison
            plot_final_comparison_variants(effective_diffusivity_thesis, theoretical_variants, ...
                constants, plots_dir);
        end
        
        if isfield(plot_flags, 'collisionality_comparison') && plot_flags.collisionality_comparison
             plot_collisionality_comparison(derived_profiles, constants, velocity_results, ...
                effective_diffusivity_thesis, plots_dir);
        end

        fprintf('--- Plotting complete. Figures saved to %s ---\n', plots_dir);
    catch ME
        warning('TCABR:PlottingError', 'Plotting phase failed: %s', ME.message);
        fprintf(2, 'Error occurred in file %s at line %d.\n', ME.stack(1).file, ME.stack(1).line);
    end
end

fprintf('\n>>> TCABR modular momentum analysis completed successfully.\n');


% =========================================================================
% LOCAL STAGE WRAPPERS
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

function theoretical_models = stage4_models(temperature_results, magnetic_field, derived_profiles, constants)
    theoretical_models = compute_theoretical_models(temperature_results, magnetic_field, derived_profiles, constants);
end

function [neutral_profile, collision_profiles] = stage5_neutral_collision(r_fine, temperature_results, constants)
    neutral_profile = compute_neutral_density_profile(r_fine, constants);
    collision_profiles = compute_collision_profiles(temperature_results, neutral_profile, r_fine);
end

function effective_diffusivity_thesis = stage6_diffusivity_thesis(velocity_results, collision_profiles, r_fine)
    effective_diffusivity_thesis = compute_diffusivity_profile(velocity_results, collision_profiles, r_fine);
end

function scaling_law_results = stage7a_diffusivity_scaling(velocity_results, derived_profiles, constants, r_fine)
    scaling_law_results = compute_chi_eff_from_scalings(velocity_results, derived_profiles, constants, r_fine);
end

% NOTE: build_theoretical_variants wrapper removed to avoid infinite recursion.
% The main script calls the function from 'src/' directly.