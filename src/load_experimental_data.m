function [exp_data] = load_experimental_data(constants)
%LOAD_EXPERIMENTAL_DATA Loads and structures the experimental data from TCABR.
%   This function performs two main tasks:
%   1. Reads the base profile data (ne, ni, etc.) from a text file. It
%      robustly handles file I/O and validates the loaded content.
%   2. Defines the specific, hard-coded experimental data points for
%      toroidal velocity and ion temperature used in subsequent analyses.
%
%   CRITICAL OPERATION: It assumes the first column in the raw data file is a
%   NORMALISED radius and converts it to physical units (meters) by multiplying
%   by the minor radius 'a'.
%
%   Syntax:
%       exp_data = load_experimental_data(constants)
%
%   Input:
%       constants - A structure containing file paths (constants.data.path)
%                   and machine parameters (constants.machine.a).
%
%   Output:
%       exp_data - A structure containing all raw and experimental data, with
%                  all vector fields ensured to be column vectors.

fprintf('Loading experimental data...\n');

a = constants.machine.a; % Minor radius for de-normalisation

% --- 1. Load Base Profiles from File ---
data_filepath = constants.data.path;
try
    if ~exist(data_filepath, 'file')
        error('Data file not found at: %s', data_filepath);
    end
    opts = detectImportOptions(data_filepath, 'FileType', 'text');
    opts.VariableNames = {'r_norm_raw', 'Te', 'Ti_unused', 'ne', 'ni', 'nCI', ...
        'nCII', 'nCIII', 'nCIV', 'nCV', 'nCVI', 'nCVII', 'n0', 'pi', 'p'};
    opts.DataLines = 2;
    raw_table = readtable(data_filepath, opts);
catch ME
    fprintf(2, 'Failed to read the experimental profile data file.\n');
    rethrow(ME);
end

% Validate loaded data
validateattributes(raw_table.r_norm_raw, {'numeric'}, {'vector', 'nonempty', 'real'}, mfilename, 'r_norm_raw');
validateattributes(raw_table.ne, {'numeric'}, {'vector', 'nonempty', 'real'}, mfilename, 'ne');
validateattributes(raw_table.ni, {'numeric'}, {'vector', 'nonempty', 'real'}, mfilename, 'ni');

% Convert normalised radius to physical units [meters] and ensure column vectors
exp_data.r_profiles = raw_table.r_norm_raw(:) * a;
exp_data.ne_raw = raw_table.ne(:); % Raw electron density profile [m^-3]
exp_data.ni_raw = raw_table.ni(:); % Raw ion density profile [m^-3]

fprintf('Base profiles loaded and radial grid converted to meters.\n');

% --- 2. Define Hard-coded Experimental Data Points ---
% These are the specific, discrete data points used for the primary fits.
% Source: [Specify source, e.g., TCABR shot #XXXXX, or publication]

% Experimental Toroidal Velocity Data [m] and [km/s]
exp_data.r_Vphi_exp = [5.02E-04; 0.0104; 0.02036; 0.03026; 0.04009; 0.05005; 0.06001; ...
    0.06991; 0.07981; 0.08977; 0.10013; 0.11; 0.11989; 0.12985; 0.13966; ...
    0.14962; 0.15958; 0.1695; 0.17971];

exp_data.Vphi_exp = -[21.03307; 24.12041; 25.06326; 27.04137; 25.96912; 24.15739; ...
    20.05325; 18.62974; 13.1021; 12.529; 10.99457; 8.12907; 8.07361; ...
    8.0921; 2.9342; -4.0909; 1.99135; 2.89722; 1.10397];

exp_data.Vphi_err_exp = abs(-[0.94284; 1.81174; 0.90586; 2.01510; 1.44200; 2.34786; ...
    1.81174; 1.42351; 1.79325; 1.42351; 0.53613; 1.23864; 0.88738; ...
    0.88738; 3.08734; 1.99660; 4.58480; 3.49406; 3.38314]);

% Experimental Ion Temperature Data [m] and [eV]
exp_data.r_Ti_exp = [0.1021; 4.5631; 6.6125; 10; 12; 14; 16; 17] / 100;
exp_data.Ti_exp_raw = [154.3163; 142.5772; 118.2958; 81; 65; 51; 31; 32];
exp_data.Ti_err_exp_raw = [12.5423; 8.897; 9.7002; 16; 15; 13; 10; 11];

% "Reconstructed" Ti data used for the canonical fit, as per the original analysis.
exp_data.Ti_reconstructed_exp = [269.9322956; 212.6290401; 154.5794833; ...
    94.59042348; 72.86028761; 58.16568356; 30.5396089; 32];
relative_errors = exp_data.Ti_err_exp_raw ./ exp_data.Ti_exp_raw;
exp_data.Ti_reconstructed_err_exp = relative_errors .* exp_data.Ti_reconstructed_exp;

% Final check to ensure all vector outputs are column vectors
fields = fieldnames(exp_data);
for i = 1:numel(fields)
    if isvector(exp_data.(fields{i})) && ~ischar(exp_data.(fields{i}))
        exp_data.(fields{i}) = exp_data.(fields{i})(:);
    end
end

fprintf('Specific experimental data points defined.\n');
fprintf('Data loading complete.\n\n');

end