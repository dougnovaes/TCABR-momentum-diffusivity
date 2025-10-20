% =========================================================================
% TCABR MOMENTUM ANALYSIS - MAIN SCRIPT (MODULAR & OPTIMISED)
% =========================================================================
% Modular, cache-based orchestrator with per-stage plotting and robust I/O.
% =========================================================================

clc; clear; close all;
fprintf('Initiating TCABR modular momentum analysis script...\n\n');

addpath('src', 'plotting', 'utils');

% Results directory
results_dir = 'results';
if ~exist(results_dir, 'dir'), mkdir(results_dir); end

%% ---------------------- CONTROL PANEL -----------------------------------
t = true; f = false;

% run_flags: true = force recompute stage; false = load if present
run_flags = struct( ...
    'setup_and_data',        t, ...
    'profile_analyses',      t, ...
    'physics_profiles',      t, ...
    'theoretical_models',    t, ...
    'neutral_and_collision', t, ...
    'diffusivity_thesis',    t, ...
    'diffusivity_scaling',   t ...
);

% plot_flags: control which stage-specific plots are shown
plot_flags = struct( ...
    'stage1_setup',        t, ...
    'stage2_profiles',     t,  ...
    'stage3_physics',      t, ...
    'stage4_models',       t,  ...
    'stage5_collisions',   t,  ...
    'stage6_diffusivity',  t,  ...
    'stage7_scalings',     t ...
);

% Behaviour: stop on critical error?
stop_on_error = true;

%% ---------------------- STAGE 1: Setup & Data ----------------------------
try
    [constants, exp_data, r_fine] = load_or_compute( ...
        fullfile(results_dir, 'stage1_setup_data.mat'), ...
        run_flags.setup_and_data, ...
        @() stage1_setup_and_data(), ...
        {'constants', 'exp_data', 'r_fine'} ...
    );
    if plot_flags.stage1_setup
        try plot_stage1_setup(constants, exp_data, r_fine); catch ME, warning(ME.identifier, 'Plot stage1 failed: %s', ME.message); end
    end
catch ME
    fprintf('Stage1 failed: %s\n', ME.message);
    if stop_on_error, rethrow(ME); end
end

%% ---------------------- STAGE 2: Profile Analyses ------------------------
try
    [velocity_results, temperature_results] = load_or_compute( ...
        fullfile(results_dir, 'stage2_profiles.mat'), ...
        run_flags.profile_analyses, ...
        @() stage2_profiles(exp_data, constants, r_fine), ...
        {'velocity_results', 'temperature_results'} ...
    );
    if plot_flags.stage2_profiles
        try plot_stage2_profiles(exp_data, velocity_results, temperature_results, r_fine, constants); catch ME, warning(ME.identifier, 'Plot stage2 failed: %s', ME.message); end
    end
catch ME
    fprintf('Stage2 failed: %s\n', ME.message);
    if stop_on_error, rethrow(ME); end
end

%% ---------------------- STAGE 3: Physics Profiles ------------------------
try
    [magnetic_field, derived_profiles] = load_or_compute( ...
        fullfile(results_dir, 'stage3_physics.mat'), ...
        run_flags.physics_profiles, ...
        @() stage3_physics(exp_data, temperature_results, constants, r_fine), ...
        {'magnetic_field', 'derived_profiles'} ...
    );
    if plot_flags.stage3_physics
        try plot_stage3_physics(magnetic_field, derived_profiles, constants, r_fine); catch ME, warning(ME.identifier, 'Plot stage3 failed: %s', ME.message); end
    end
catch ME
    fprintf('Stage3 failed: %s\n', ME.message);
    if stop_on_error, rethrow(ME); end
end

%% ---------------------- STAGE 4: Theoretical Models ---------------------
try
    theoretical_models = load_or_compute( ...
        fullfile(results_dir, 'stage4_models.mat'), ...
        run_flags.theoretical_models, ...
        @() stage4_models(temperature_results, magnetic_field, derived_profiles, constants), ...
        {'theoretical_models'} ...
    );
    if plot_flags.stage4_models
        try plot_stage4_models(theoretical_models, derived_profiles, constants, r_fine); catch ME, warning(ME.identifier, 'Plot stage4 failed: %s', ME.message); end
    end
catch ME
    fprintf('Stage4 failed: %s\n', ME.message);
    if stop_on_error, rethrow(ME); end
end

%% ---------------------- STAGE 5: Neutral & Collisions -------------------
try
    [neutral_profile, collision_profiles] = load_or_compute( ...
        fullfile(results_dir, 'stage5_collisions.mat'), ...
        run_flags.neutral_and_collision, ...
        @() stage5_neutral_collision(r_fine, temperature_results, constants), ...
        {'neutral_profile', 'collision_profiles'} ...
    );
    if plot_flags.stage5_collisions
        try plot_stage5_collisions(neutral_profile, collision_profiles, r_fine, constants); catch ME, warning(ME.identifier, 'Plot stage5 failed: %s', ME.message); end
    end
catch ME
    fprintf('Stage5 failed: %s\n', ME.message);
    if stop_on_error, rethrow(ME); end
end

%% ---------------------- STAGE 6: Effective Diffusivity (Thesis) --------
try
    effective_diffusivity_thesis = load_or_compute( ...
        fullfile(results_dir, 'stage6_diffusivity_thesis.mat'), ...
        run_flags.diffusivity_thesis, ...
        @() stage6_diffusivity_thesis(velocity_results, collision_profiles, r_fine), ...
        {'effective_diffusivity_thesis'} ...
    );
    if plot_flags.stage6_diffusivity
        try plot_stage6_diffusivity(effective_diffusivity_thesis, r_fine, constants); catch ME, warning(ME.identifier, 'Plot stage6 failed: %s', ME.message); end
    end
catch ME
    fprintf('Stage6 failed: %s\n', ME.message);
    if stop_on_error, rethrow(ME); end
end

%% ---------------------- STAGE 7: Effective Diffusivity (Scaling) -------
try
    scaling_law_results = load_or_compute( ...
        fullfile(results_dir, 'stage7_diffusivity_scaling.mat'), ...
        run_flags.diffusivity_scaling, ...
        @() stage7_diffusivity_scaling(velocity_results, derived_profiles, constants, r_fine), ...
        {'scaling_law_results'} ...
    );
    if plot_flags.stage7_scalings
        try plot_stage7_scalings(scaling_law_results, r_fine, constants); catch ME, warning(ME.identifier, 'Plot stage7 failed: %s', ME.message); end
    end
catch ME
    fprintf('Stage7 failed: %s\n', ME.message);
    if stop_on_error, rethrow(ME); end
end

fprintf('\nAll requested stages processed.\n\n');

%% ---------------------- FINAL: Comparison Plot --------------------------
% Always optional — call only when results exist
if exist('effective_diffusivity_thesis', 'var') && exist('scaling_law_results', 'var')
    try
        plot_diffusivity_comparison(effective_diffusivity_thesis, scaling_law_results, constants, r_fine);
    catch ME
        warning(ME.identifier, 'Final comparison plot failed: %s', ME.message);
    end
end

fprintf('Script finished.\n');

%% =======================================================================
%% ------------------ load_or_compute helper (optimised) ------------------
%% =======================================================================
function varargout = load_or_compute(filepath, force_recalc, compute_func, out_names)
% LOAD_OR_COMPUTE  Load cached results or compute and save them.
% Usage:
%   [a,b,...] = load_or_compute(filepath, force_recalc, @() compute(), {'a','b',...})
%
% If out_names is provided, the saved .mat will contain those variable names.
% The function returns as many outputs as the caller expects.

if nargin < 2 || isempty(force_recalc), force_recalc = false; end
if nargin < 3, error('load_or_compute requires (filepath, force_recalc, compute_func).'); end
if nargin < 4, out_names = {}; end

% Normalise filepath
if isstring(filepath), filepath = char(filepath); end
[dirpath, fname, fext] = fileparts(filepath);
if isempty(fext), filepath = [filepath '.mat']; fext = '.mat'; end
if ~isempty(dirpath) && ~exist(dirpath, 'dir'), mkdir(dirpath); end

% Decide load or compute
if exist(filepath, 'file') && ~force_recalc
    S = load(filepath);
    fields = fieldnames(S);
    nout = nargout;
    % If out_names provided and exist in file, return in that order
    if ~isempty(out_names) && all(ismember(out_names, fields))
        for k = 1:nout
            varargout{k} = S.(out_names{k});
        end
    else
        % fallback: try canonical 'out1..outN'
        canonical = arrayfun(@(k) sprintf('out%d', k), 1:nout, 'UniformOutput', false);
        if all(ismember(canonical, fields))
            for k = 1:nout
                varargout{k} = S.(canonical{k});
            end
        else
            % final fallback: return first nout fields (warn)
            if numel(fields) < nout
                error('File %s contains %d fields but %d outputs requested.', filepath, numel(fields), nout);
            end
            for k = 1:nout
                varargout{k} = S.(fields{k});
            end
            warning('Loaded %s: returning first %d fields (%s). Consider re-saving with explicit out_names.', filepath, nout, strjoin(fields(1:nout), ', '));
        end
    end
    fprintf('Loaded cached results: %s\n', filepath);
    return;
end

% Compute
fprintf('Computing stage: %s\n', [fname fext]);
[varargout{1:nargout}] = compute_func();

% Prepare struct for saving
S = struct();
if ~isempty(out_names) && numel(out_names) >= nargout
    for k = 1:nargout
        name = out_names{k};
        S.(name) = varargout{k};
    end
else
    for k = 1:nargout
        S.(sprintf('out%d', k)) = varargout{k};
    end
end

% Save
try
    save(filepath, '-struct', 'S');
    fprintf('Saved results: %s\n', filepath);
catch ME
    warning('Failed to save %s: %s', filepath, ME.message);
end
end

%% =======================================================================
%% ------------------------ Stage wrappers (unchanged) ---------------------
%% =======================================================================
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
