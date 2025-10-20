% =========================================================================
% TCABR MOMENTUM ANALYSIS - MAIN SCRIPT (MODULAR VERSION)
% =========================================================================
% Modular and cache-based orchestrator for the analysis of toroidal
% momentum transport in the TCABR tokamak.
% Each stage saves and loads its own results independently.
% =========================================================================

clc; clear; close all;
fprintf('Initiating TCABR modular momentum analysis script...\n\n');

addpath('src', 'plotting', 'utils');
results_dir = 'results';
if ~exist(results_dir, 'dir'), mkdir(results_dir); end

%% --- EXECUTION CONTROL PANEL --------------------------------------------
% Define which stages to run or reload from cache.
% true  -> force recalculation of that stage
% false -> load from saved .mat if available
run_flags = struct( ...
    'setup_and_data',           false, ...
    'profile_analyses',         true, ...
    'physics_profiles',         true, ...
    'theoretical_models',       true, ...
    'neutral_and_collision',    true, ...
    'diffusivity_thesis',       true, ...
    'diffusivity_scaling',      true ...
);

%% --- STAGE 1: Setup and Data Loading -----------------------------------
[constants, exp_data, r_fine] = load_or_compute( ...
    fullfile(results_dir, 'stage1_setup_data.mat'), ...
    run_flags.setup_and_data, ...
    @() stage1_setup_and_data() ...
);

%% --- STAGE 2: Profile Analyses -----------------------------------------
[velocity_results, temperature_results] = load_or_compute( ...
    fullfile(results_dir, 'stage2_profiles.mat'), ...
    run_flags.profile_analyses, ...
    @() stage2_profiles(exp_data, constants, r_fine) ...
);

%% --- STAGE 3: Physics Profile Calculations ------------------------------
[magnetic_field, derived_profiles] = load_or_compute( ...
    fullfile(results_dir, 'stage3_physics.mat'), ...
    run_flags.physics_profiles, ...
    @() stage3_physics(exp_data, temperature_results, constants, r_fine) ...
);

%% --- STAGE 4: Theoretical Model Calculations ----------------------------
theoretical_models = load_or_compute( ...
    fullfile(results_dir, 'stage4_models.mat'), ...
    run_flags.theoretical_models, ...
    @() stage4_models(temperature_results, magnetic_field, derived_profiles, constants) ...
);

%% --- STAGE 5: Neutral Density and Collisions ----------------------------
[neutral_profile, collision_profiles] = load_or_compute( ...
    fullfile(results_dir, 'stage5_collisions.mat'), ...
    run_flags.neutral_and_collision, ...
    @() stage5_neutral_collision(r_fine, temperature_results, constants) ...
);

%% --- STAGE 6: Effective Diffusivity (Thesis Method) --------------------
effective_diffusivity_thesis = load_or_compute( ...
    fullfile(results_dir, 'stage6_diffusivity_thesis.mat'), ...
    run_flags.diffusivity_thesis, ...
    @() stage6_diffusivity_thesis(velocity_results, collision_profiles, r_fine) ...
);

%% --- STAGE 7: Effective Diffusivity (Scaling Laws) ---------------------
scaling_law_results = load_or_compute( ...
    fullfile(results_dir, 'stage7_diffusivity_scaling.mat'), ...
    run_flags.diffusivity_scaling, ...
    @() stage7_diffusivity_scaling(velocity_results, derived_profiles, constants, r_fine) ...
);

fprintf('\nAll calculations complete.\n\n');

%% --- STAGE 8: Plotting -------------------------------------------------
plot_diffusivity_comparison(effective_diffusivity_thesis, scaling_law_results, constants, r_fine);
% run_verification_plots(...);

fprintf('Script finished.\n');



% ========================================================================
% ======================== AUXILIARY FUNCTIONS ===========================
% ========================================================================

function varargout = load_or_compute(filepath, force_recalc, compute_func)
% LOAD_OR_COMPUTE  Load cached results or compute-and-save.
% Usage:
%   [a,b,...] = load_or_compute(filepath, force_recalc, @() compute_stage());
%
% - filepath: full path to .mat file to load/save
% - force_recalc: logical; if true, compute even if file exists
% - compute_func: function handle that returns the outputs for this stage
%
% This function always returns as many outputs as the caller expects.

if nargin < 2 || isempty(force_recalc)
    force_recalc = false;
end
if nargin < 3
    error('load_or_compute requires (filepath, force_recalc, compute_func).');
end

% Normalise filepath to char
if isstring(filepath), filepath = char(filepath); end
[~, fname, fext] = fileparts(filepath);
if isempty(fext)
    filepath = [filepath '.mat'];
end

% Decide: load or compute
if exist(filepath, 'file') && ~force_recalc
    try
        S = load(filepath);
    catch ME
        error('Failed to load %s: %s', filepath, ME.message);
    end

    % If file contains out1...outN (our canonical save), return them in order.
    outN = nargout;
    outFields = fieldnames(S);
    hasCanonical = all(arrayfun(@(k) ismember(sprintf('out%d', k), outFields), 1:outN));

    if hasCanonical
        for k = 1:outN
            varargout{k} = S.(sprintf('out%d', k));
        end
    else
        % Fallback: if number of fields equals requested outputs, return in field order
        if numel(outFields) >= outN
            for k = 1:outN
                varargout{k} = S.(outFields{k});
            end
            warning('Loaded %s: fields returned in stored order (%s).', filepath, strjoin(outFields(1:outN), ', '));
        else
            error('File %s does not contain enough outputs (%d requested, %d available).', filepath, outN, numel(outFields));
        end
    end

    fprintf('Loaded cached results: %s\n', filepath);
    return;
end

% --- compute and save ---
fprintf('Computing stage: %s\n', [fname '.mat']);
[varargout{1:nargout}] = compute_func();

% Prepare structure to save: canonical names out1,out2,...
S = struct();
for k = 1:nargout
    S.(sprintf('out%d', k)) = varargout{k};
end

% Create directory if needed
[outdir, ~, ~] = fileparts(filepath);
if ~isempty(outdir) && ~exist(outdir, 'dir')
    mkdir(outdir);
end

try
    save(filepath, '-struct', 'S');
    fprintf('Saved results: %s\n', filepath);
catch ME
    warning('Failed to save results to %s: %s', filepath, ME.message);
end
end

function names = get_output_names(n)
    % Default output names to avoid empty fieldnames
    names = arrayfun(@(k) sprintf('out%d', k), 1:n, 'UniformOutput', false);
end


% ========================================================================
% ===================== STAGE FUNCTION WRAPPERS ==========================
% ========================================================================

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