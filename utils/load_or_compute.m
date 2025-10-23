function varargout = load_or_compute(filepath, force_recalc, compute_func, out_names)
%LOAD_OR_COMPUTE Loads a .mat file or executes a function, with robust cache validation.
%   This function manages the caching mechanism for the analysis pipeline.
%
%   BEHAVIOR:
%   1. Checks if 'filepath' exists and 'force_recalc' is false.
%   2. If so, it loads the .mat file and validates that ALL variable names
%      specified in 'out_names' exist as fields in the loaded data.
%   3. If the cache is valid, it returns the loaded variables.
%   4. If the file does not exist, 'force_recalc' is true, or the cache is
%      invalid (stale), it executes the 'compute_func' handle.
%   5. The results of 'compute_func' are then saved to 'filepath' and returned.
%
%   Syntax:
%       [a, b] = load_or_compute('results/file.mat', false, @my_func, {'a', 'b'})
%
%   Inputs:
%       filepath     - Full path to the target .mat cache file.
%       force_recalc - Logical. If true, ignores existing cache and forces re-computation.
%       compute_func - A function handle with no arguments that returns the variables to be cached.
%       out_names    - Cell array of strings with the names of the output variables.
%                      The order must match the output order of 'compute_func'.
%
%   Outputs:
%       varargout    - The results, either from the cache or from computation.

    arguments
        filepath (1,1) string
        force_recalc (1,1) logical
        compute_func (1,1) function_handle
        out_names (1,:) cell
    end

    % --- 1. Check for valid cache ---
    if exist(filepath, 'file') && ~force_recalc
        fprintf('Found cache: %s. Validating...\n', filepath);
        S = load(filepath);
        
        % Validate that all expected fields exist in the loaded struct
        if all(isfield(S, out_names))
            fprintf('...Cache is valid. Loading results.\n');
            varargout = cell(1, nargout);
            for k = 1:nargout
                varargout{k} = S.(out_names{k});
            end
            return;
        else
            fprintf('...Cache is stale (missing variables). Re-computing.\n');
        end
    end

    % --- 2. Compute, as cache is absent, invalid, or ignored ---
    fprintf('Computing and saving results to: %s\n', filepath);

    % Ensure the target directory exists before trying to save
    [target_dir, ~, ~] = fileparts(filepath);
    if ~isempty(target_dir) && ~exist(target_dir, 'dir')
        mkdir(target_dir);
    end

    % Execute the computation function
    [varargout{1:nargout}] = compute_func();

    % --- 3. Save the new results to cache ---
    S = struct();
    for k = 1:nargout
        S.(out_names{k}) = varargout{k};
    end
    
    try
        save(filepath, '-struct', 'S', '-v7.3');
        fprintf('...Saved results successfully.\n');
    catch ME
        warning('Could not save cache file: %s\nMessage: %s', filepath, ME.message);
    end
end