% =========================================================================
% main_analysis.m
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
%     ├── src/         # Core computational routines (physics & analysis)
%     ├── utils/       # Utility functions and wrappers
%     ├── plotting/    # Publication-quality visualisation scripts
%     ├── studies/     # High-level parameter scans or specialised analyses
%     └── results/     # Per-stage .mat results and plots
%
% REQUIREMENTS:
%   MATLAB R2025b or newer
% =========================================================================

clc; clear; close all;
fprintf('>>> Initiating TCABR modular momentum analysis pipeline...\n\n');

% -------------------------------------------------------------------------
% PATH SETUP
% -------------------------------------------------------------------------
addpath('src', 'plotting', 'utils');

results_dir = 'results';
if ~exist(results_dir, 'dir'), mkdir(results_dir); end
plots_dir = fullfile(results_dir, 'plots');
if ~exist(plots_dir, 'dir'), mkdir(plots_dir); end

% -------------------------------------------------------------------------
% CONTROL PANEL
% -------------------------------------------------------------------------
t = true; f = false;  % convenience flags

run_flags = struct( ...
    'setup_and_data',        f, ...
    'profile_analyses',      f, ...
    'physics_profiles',      f, ...
    'theoretical_models',    f, ...
    'neutral_and_collision', f, ...
    'diffusivity_thesis',    f, ...
    'diffusivity_scaling',   f ...
);

plot_flags = struct( ...
    'justification_and_fits',   t, ...
    'thesis_method_results',    f, ...
    'supporting_profiles',      f, ...
    'final_comparison',         f ...
);

study_flags = struct( ...
    'neutral_scan', f ... % activate studies/scan_neutrals_fit
);

stop_on_error = true;  % if true, halts execution upon any unhandled error

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

%% -------------------- STAGE 7: DIFFUSIVITY (SCALING LAWS) ---------------
% Outputs: scaling_law_results
try
    scaling_law_results = load_or_compute( ...
        fullfile(results_dir, 'stage7_diffusivity_scaling.mat'), ...
        run_flags.diffusivity_scaling, ...
        @() stage7_diffusivity_scaling(velocity_results, derived_profiles, constants, r_fine), ...
        {'scaling_law_results'});
    fprintf('[✓] Stage 7 completed: effective diffusivity (scaling laws).\n');
catch ME
    handle_stage_error(ME, 7, stop_on_error);
end

fprintf('\n>>> All calculation stages executed (or loaded) successfully.\n\n');

% =========================================================================
% PLOTTING SECTION
% =========================================================================
fprintf('--- Generating selected plots ---\n');

try
    if plot_flags.justification_and_fits
        plot_justification_and_fits(constants, exp_data, velocity_results, ...
            temperature_results, theoretical_models, r_fine);
    end

    if plot_flags.thesis_method_results
        plot_thesis_method_results(constants, neutral_profile, collision_profiles, ...
            effective_diffusivity_thesis, r_fine);
    end

    if plot_flags.supporting_profiles
        plot_supporting_profiles(constants, derived_profiles, r_fine);
    end

    if plot_flags.final_comparison
        % Compute or load chi_eff variants for theoretical comparison
        if ~exist('variants','var') || isempty(variants)
            variants = compute_chi_eff_variants(scaling_law_results, ...
                velocity_results, derived_profiles, constants, r_fine);
        end
        plot_final_comparison(effective_diffusivity_thesis, variants, ...
            scaling_law_results, constants);
    end
catch ME
    warning(ME.identifier,'Plotting phase failed: %s', ME.message);
end

% =========================================================================
% STUDY SECTION (optional)
% =========================================================================
if study_flags.neutral_scan
    try
        fprintf('--- Running neutral density parameter scan ---\n');
        addpath('studies');
        run('studies/run_scan_neutrals_fit.m');
    catch ME
        warning(ME.identifier,'Neutral scan study failed: %s', ME.message);
    end
end

fprintf('\n>>> TCABR modular momentum analysis completed successfully.\n');

% =========================================================================
% LOCAL UTILITY FUNCTIONS
% =========================================================================
function handle_stage_error(ME, stage_num, stop_on_error)
    fprintf(2, '[X] Stage %d failed: %s\n', stage_num, ME.message);
    if stop_on_error, rethrow(ME); end
end

function varargout = load_or_compute(filepath, force_recalc, compute_func, out_names)
    if exist(filepath, 'file') && ~force_recalc
        fprintf('Loading cached results from: %s\n', filepath);
        S = load(filepath);
        varargout = cell(1, nargout);
        for k = 1:nargout
            varargout{k} = S.(out_names{k});
        end
        return;
    end
    fprintf('Computing and saving results to: %s\n', filepath);
    [varargout{1:nargout}] = compute_func();
    S = struct();
    for k = 1:nargout
        S.(out_names{k}) = varargout{k};
    end
    save(filepath, '-struct', 'S');
end

% --- Stage wrappers (simple delegation to /src functions) ----------------
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

function scaling_law_results = stage7_diffusivity_scaling(velocity_results, derived_profiles, constants, r_fine)
    scaling_law_results = compute_chi_eff_from_scalings(velocity_results, derived_profiles, constants, r_fine);
end

% =========================================================================
% End of main_analysis.m
% =========================================================================
