% =========================================================================
% TCABR MOMENTUM TRANSPORT ANALYSIS (MODULAR & CACHE-ENABLED)
% =========================================================================
% PURPOSE:
%   This main script orchestrates the entire data analysis pipeline for
%   toroidal momentum transport in the TCABR tokamak. It is designed for
%   modular execution, per-stage caching, and controlled plotting.
%
% DESIGN:
%   - Modular "stage" system: results saved to /results as .mat per stage.
%   - Clear separation between computation, visualisation, and studies.
%   - Interruptible pipeline via 'run_flags' and 'plot_flags'.
%
% LAYOUT:
%   Root/
%     ├── main_analysis.m    # This script
%     ├── data/              # Input data files (e.g., .txt)
%     ├── src/               # Core computational routines (physics & analysis)
%     ├── utils/             # Utility functions and wrappers
%     ├── plotting/          # Publication-quality visualisation scripts
%     ├── studies/           # High-level parameter scans or specialised analyses
%     └── results/           # Per-stage .mat results and generated plots
%
% REQUIREMENTS:
%   MATLAB R2025b or newer
%   - Curve Fitting Toolbox
%   - Parallel Computing Toolbox
%   - Statistics and Machine Learning Toolbox
% =========================================================================

clc; clear; close all;
fprintf('>>> Initiating TCABR modular momentum analysis pipeline...\n\n');

% -------------------------------------------------------------------------
% PATH SETUP (Robust)
% -------------------------------------------------------------------------
[base_dir, ~, ~] = fileparts(mfilename('fullpath'));
addpath(fullfile(base_dir, 'src'), ...
        fullfile(base_dir, 'plotting'), ...
        fullfile(base_dir, 'utils'), ...
        fullfile(base_dir, 'studies'));

results_dir = fullfile(base_dir, 'results');
if ~exist(results_dir, 'dir'), mkdir(results_dir); end
plots_dir = fullfile(results_dir, 'plots');
if ~exist(plots_dir, 'dir'), mkdir(plots_dir); end

% -------------------------------------------------------------------------
% CONTROL PANEL
% -------------------------------------------------------------------------
% Set flags to 'true' to force re-computation or generate plots.
% Set to 'false' to load from cache or skip plotting.

run_flags = struct( ...
    'setup_and_data',        false, ... % Stage 1
    'profile_analyses',      false, ... % Stage 2
    'physics_profiles',      true, ... % Stage 3
    'theoretical_models',    true, ... % Stage 4
    'neutral_and_collision', false, ... % Stage 5
    'diffusivity_thesis',    false, ... % Stage 6
    'diffusivity_scaling',   true, ... % Stage 7a
    'build_variants',        true  ... % Stage 7b
);

plot_flags = struct( ...
    'justification_and_fits',  true, ...
    'thesis_method_results',   true, ...
    'supporting_profiles',     true, ...
    'final_comparison',        true,  ... % Will use the new variants
    'collisionality_comparison', true  ... % <-- NEW FLAG ADDED
);

study_flags = struct( ...
    'neutral_scan', false ... % Activates studies/run_scan_neutrals_fit
);

stop_on_error = true; % if true, halts execution upon any unhandled error

fprintf('>>> Control panel configured.\n\n');

% =========================================================================
% PIPELINE EXECUTION
% =========================================================================
% Each stage below either loads existing cached results or recomputes and
% saves them to the results directory.

%% -------------------- STAGE 1: SETUP & DATA ------------------------------
% Outputs: constants, exp_data, r_fine
try
    [constants, exp_data, r_fine] = load_or_compute( ...
        fullfile(results_dir, 'stage1_setup_data.mat'), ...
        run_flags.setup_and_data, ...
        @stage1_setup_and_data, ...
        {'constants','exp_data','r_fine'});
    fprintf('[✓] Stage 1 completed: setup & data.\n');
catch ME
    handle_stage_error(ME, 1, stop_on_error);
end

%% -------------------- STAGE 2: PROFILE ANALYSES --------------------------
% Outputs: velocity_results, temperature_results
try
    [velocity_results, temperature_results] = load_or_compute( ...
        fullfile(results_dir, 'stage2_profiles.mat'), ...
        run_flags.profile_analyses, ...
        @() stage2_profiles(exp_data, constants, r_fine), ...
        {'velocity_results','temperature_results'});
    fprintf('[✓] Stage 2 completed: profile analyses.\n');
catch ME
    handle_stage_error(ME, 2, stop_on_error);
end

%% -------------------- STAGE 3: PHYSICS PROFILES --------------------------
% Outputs: magnetic_field, derived_profiles
try
    [magnetic_field, derived_profiles] = load_or_compute( ...
        fullfile(results_dir, 'stage3_physics.mat'), ...
        run_flags.physics_profiles, ...
        @() stage3_physics(exp_data, temperature_results, constants, r_fine), ...
        {'magnetic_field','derived_profiles'});
    fprintf('[✓] Stage 3 completed: physics profiles.\n');
catch ME
    handle_stage_error(ME, 3, stop_on_error);
end

%% -------------------- STAGE 4: THEORETICAL MODELS -----------------------
% Outputs: theoretical_models
try
    theoretical_models = load_or_compute( ...
        fullfile(results_dir, 'stage4_models.mat'), ...
        run_flags.theoretical_models, ...
        @() stage4_models(temperature_results, magnetic_field, derived_profiles, constants), ...
        {'theoretical_models'});
    fprintf('[✓] Stage 4 completed: theoretical models.\n');
catch ME
    handle_stage_error(ME, 4, stop_on_error);
end

%% -------------------- STAGE 5: NEUTRALS & COLLISIONS --------------------
% Outputs: neutral_profile, collision_profiles
try
    [neutral_profile, collision_profiles] = load_or_compute( ...
        fullfile(results_dir, 'stage5_collisions.mat'), ...
        run_flags.neutral_and_collision, ...
        @() stage5_neutral_collision(r_fine, temperature_results, constants), ...
        {'neutral_profile','collision_profiles'});
    fprintf('[✓] Stage 5 completed: neutrals & collisions.\n');
catch ME
    handle_stage_error(ME, 5, stop_on_error);
end

%% -------------------- STAGE 6: DIFFUSIVITY (THESIS METHOD) --------------
% Outputs: effective_diffusivity_thesis
try
    effective_diffusivity_thesis = load_or_compute( ...
        fullfile(results_dir, 'stage6_diffusivity_thesis.mat'), ...
        run_flags.diffusivity_thesis, ...
        @() stage6_diffusivity_thesis(velocity_results, collision_profiles, r_fine), ...
        {'effective_diffusivity_thesis'});
    fprintf('[✓] Stage 6 completed: effective diffusivity (thesis method).\n');
catch ME
    handle_stage_error(ME, 6, stop_on_error);
end

%% -------------------- STAGE 7A: DIFFUSIVITY (SCALING LAWS - BASE) --------
% Outputs: scaling_law_results
try
    scaling_law_results = load_or_compute( ...
        fullfile(results_dir, 'stage7a_diffusivity_scaling_base.mat'), ...
        run_flags.diffusivity_scaling, ...
        @() stage7a_diffusivity_scaling(velocity_results, derived_profiles, constants, r_fine), ...
        {'scaling_law_results'});
    fprintf('[✓] Stage 7a completed: effective diffusivity (base scaling laws).\n');
catch ME
    handle_stage_error(ME, 7, stop_on_error);
end

%% -------------------- STAGE 7B: BUILD THEORETICAL VARIANTS ---------------
% Outputs: theoretical_variants
try
    theoretical_variants = load_or_compute( ...
        fullfile(results_dir, 'stage7b_theoretical_variants.mat'), ...
        run_flags.build_variants, ...
        @() build_theoretical_variants(velocity_results, ... % Chamada corrigida
                                       derived_profiles, constants, r_fine), ...
        {'theoretical_variants'});
    fprintf('[✓] Stage 7b completed: theoretical chi_eff variants built.\n');
catch ME
    handle_stage_error(ME, 7.5, stop_on_error);
end

fprintf('\n>>> All calculation stages executed (or loaded) successfully.\n\n');

% =========================================================================
% PLOTTING SECTION
% =========================================================================
plot_flags_as_cell = struct2cell(plot_flags);   % Step 1: Store into a variable
if any([plot_flags_as_cell{:}])                 % Step 2: Using the variable
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
            % This new plot function uses the detailed variants
            plot_final_comparison_variants(effective_diffusivity_thesis, theoretical_variants, ...
                constants, plots_dir);
        end

        if plot_flags.collisionality_comparison
            % Esta função precisa dos 'derived_profiles' e dos 'velocity_results'
            % e do resultado experimental 'effective_diffusivity_thesis'
            plot_collisionality_comparison(derived_profiles, constants, velocity_results, ...
                effective_diffusivity_thesis, plots_dir);
        end
        fprintf('--- Plotting complete. Figures saved to %s ---\n', plots_dir);
    catch ME
        warning(ME.identifier,'Plotting phase failed: %s', ME.message);
        fprintf(2, 'Error occurred in file %s at line %d.\n', ME.stack(1).file, ME.stack(1).line);
    end
end

% =========================================================================
% STUDY SECTION (optional)
% =========================================================================
study_flags_as_cell = struct2cell(study_flags);
if any([study_flags_as_cell{:}])
    fprintf('\n--- Running selected studies ---\n');
    if study_flags.neutral_scan
        try
            run_scan_neutrals_fit(velocity_results, temperature_results, theoretical_variants, ...
                constants, r_fine, results_dir);
        catch ME
            warning(ME.identifier, 'Neutral scan study failed: %s', ME.message);
            fprintf(2, 'Error occurred in file %s at line %d.\n', ME.stack(1).file, ME.stack(1).line);
        end
    end
end

fprintf('\n>>> TCABR modular momentum analysis completed successfully.\n');


% =========================================================================
% LOCAL STAGE WRAPPERS
% These are thin wrappers that delegate the actual work to functions in /src
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
    % Note: Renamed to clarify this is the 'base' calculation
    scaling_law_results = compute_chi_eff_from_scalings_base(velocity_results, derived_profiles, constants, r_fine);
end