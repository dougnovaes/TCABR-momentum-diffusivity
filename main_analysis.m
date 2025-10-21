% =========================================================================
% TCABR MOMENTUM ANALYSIS - MAIN SCRIPT (MODULAR & CACHE-ENABLED)
% =========================================================================
% This script orchestrates the full analysis pipeline for toroidal momentum
% transport in the TCABR tokamak. It features a per-stage caching system
% to avoid re-computation and a control panel for selective plotting.
% =========================================================================

clc; clear; close all;
fprintf('Initiating TCABR modular momentum analysis script...\n\n');

addpath('src', 'plotting', 'utils');

results_dir = 'results';
if ~exist(results_dir, 'dir'), mkdir(results_dir); end

%% ---------------------- CONTROL PANEL -----------------------------------
t = true; f = false;

run_flags = struct( ...
    'setup_and_data',        f, ...
    'profile_analyses',      f, ...
    'physics_profiles',      f, ... % << Ensure this is false unless needed
    'theoretical_models',    f, ...
    'neutral_and_collision', f, ...
    'diffusivity_thesis',    f, ...
    'diffusivity_scaling',   f ...
);

% --- Plotting Flags (Updated Names) ---
plot_flags = struct( ...
    'justification_and_fits', f, ...
    'thesis_method_results',  f, ...
    'supporting_profiles',    t, ... % << Set TRUE to run this group
    'final_comparison',       f ...
);

stop_on_error = true;

%% ========================================================================
%  --- EXECUTION PIPELINE (Calculations) ---
% =========================================================================
% (As seções de cálculo permanecem as mesmas das versões anteriores)

%% --- STAGE 1: Setup & Data ---
try
    [constants, exp_data, r_fine] = load_or_compute( ...
        fullfile(results_dir, 'stage1_setup_data.mat'), run_flags.setup_and_data, ...
        @stage1_setup_and_data, {'constants', 'exp_data', 'r_fine'});
catch ME, if stop_on_error, rethrow(ME); else, warning(ME.message); end, end

%% --- STAGE 2: Profile Analyses ---
try
    [velocity_results, temperature_results] = load_or_compute( ...
        fullfile(results_dir, 'stage2_profiles.mat'), run_flags.profile_analyses, ...
        @() stage2_profiles(exp_data, constants, r_fine), {'velocity_results', 'temperature_results'});
catch ME, if stop_on_error, rethrow(ME); else, warning(ME.message); end, end

%% --- STAGE 3: Physics Profiles ---
try
    [magnetic_field, derived_profiles] = load_or_compute( ...
        fullfile(results_dir, 'stage3_physics.mat'), run_flags.physics_profiles, ...
        @() stage3_physics(exp_data, temperature_results, constants, r_fine), {'magnetic_field', 'derived_profiles'});
catch ME, if stop_on_error, rethrow(ME); else, warning(ME.message); end, end

%% --- STAGE 4: Theoretical Models ---
try
    theoretical_models = load_or_compute( ...
        fullfile(results_dir, 'stage4_models.mat'), run_flags.theoretical_models, ...
        @() stage4_models(temperature_results, magnetic_field, derived_profiles, constants), {'theoretical_models'});
catch ME, if stop_on_error, rethrow(ME); else, warning(ME.message); end, end

%% --- STAGE 5: Neutral & Collisions ---
try
    [neutral_profile, collision_profiles] = load_or_compute( ...
        fullfile(results_dir, 'stage5_collisions.mat'), run_flags.neutral_and_collision, ...
        @() stage5_neutral_collision(r_fine, temperature_results, constants), {'neutral_profile', 'collision_profiles'});
catch ME, if stop_on_error, rethrow(ME); else, warning(ME.message); end, end

%% --- STAGE 6: Effective Diffusivity (Thesis Method) ---
try
    effective_diffusivity_thesis = load_or_compute( ...
        fullfile(results_dir, 'stage6_diffusivity_thesis.mat'), run_flags.diffusivity_thesis, ...
        @() stage6_diffusivity_thesis(velocity_results, collision_profiles, r_fine), {'effective_diffusivity_thesis'});
catch ME, if stop_on_error, rethrow(ME); else, warning(ME.message); end, end

%% --- STAGE 7: Effective Diffusivity (from Scaling Laws) ---
try
    scaling_law_results = load_or_compute( ...
        fullfile(results_dir, 'stage7_diffusivity_scaling.mat'), run_flags.diffusivity_scaling, ...
        @() stage7_diffusivity_scaling(velocity_results, derived_profiles, constants, r_fine), {'scaling_law_results'});
catch ME, if stop_on_error, rethrow(ME); else, warning(ME.message); end, end

fprintf('\nAll calculation stages processed successfully.\n\n');

%% ========================================================================
%  --- PLOTTING STAGE ---
% =========================================================================
fprintf('--- Generating Selected Plots ---\n');

if plot_flags.justification_and_fits
    try
        plot_stage1_justification_and_fits(constants, exp_data, velocity_results, temperature_results, theoretical_models, r_fine);
    catch ME, warning(ME.identifier,'Plotting failed for "Justification & Fits": %s', ME.message); end
end

if plot_flags.thesis_method_results
    try
        plot_stage2_thesis_method(constants, neutral_profile, collision_profiles, effective_diffusivity_thesis, r_fine);
    catch ME, warning(ME.identifier,'Plotting failed for "Thesis Method Results": %s', ME.message); end
end

if plot_flags.supporting_profiles
    try
        plot_stage3_supporting_profiles(constants, derived_profiles, r_fine);
    catch ME
        warning(ME.identifier,'Plotting failed for "Supporting Profiles": %s', ME.message);
    end
end

fprintf('\nScript finished.\n');

%% ========================================================================
%  --- HELPER AND WRAPPER FUNCTIONS ---
% =========================================================================
% (As wrapper functions permanecem as mesmas)

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
    fprintf('Computing and saving results to: %s\n', [fileparts(filepath), filesep, out_names{:}]);
    [varargout{1:nargout}] = compute_func();
    S = struct();
    for k = 1:nargout
        S.(out_names{k}) = varargout{k};
    end
    save(filepath, '-struct', 'S');
end

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
    effective_diffusivity_thesis = compute_effective_diffusivity(velocity_results, collision_profiles, r_fine);
end
function scaling_law_results = stage7_diffusivity_scaling(velocity_results, derived_profiles, constants, r_fine)
    scaling_law_results = compute_chi_eff_from_scalings(velocity_results, derived_profiles, constants, r_fine);
end