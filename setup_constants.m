function constants = setup_constants()
%SETUP_CONSTANTS Defines and returns a structure with all project constants.
%   This function organises all physical constants, machine parameters, analysis
%   settings, and model parameters into a nested structure for clarity and
%   easy access throughout the analysis pipeline.

% --- Machine and Discharge Parameters (TCABR) ---
constants.machine.a = 0.18;         % Minor radius [m]
constants.machine.R0 = 0.61;        % Major radius [m]
constants.machine.Ip = 90E3;        % Plasma current [A]
constants.machine.B0 = 1.07;        % On-axis magnetic field [T]
constants.machine.epsilon_aspect_ratio = constants.machine.a / constants.machine.R0; % Inverse aspect ratio

% --- Fundamental Physical Constants (SI units) ---
constants.physics.e_charge = 1.60217663E-19;      % Elementary charge [C]
constants.physics.me = 9.1093837E-31;             % Electron mass [kg]
constants.physics.mp = 1.672622E-27;             % Proton mass [kg]
constants.physics.mu0 = 4 * pi * 1E-7;            % Vacuum permeability [H/m]
constants.physics.epsilon0 = 8.854187812E-12;     % Vacuum permittivity [F/m]

% --- Plasma Species and Analysis Parameters ---
% NOTE: 'ion' refers to the main hydrogen species (H+).
% Toroidal velocity is measured from C6+ impurities, assuming V_phi_C6+ ~ V_phi_H+.
constants.plasma.m_ion = constants.physics.mp;      % Mass of main ion (Hydrogen) [kg]
constants.plasma.m_ion_amu = 1.00784;             % Atomic mass of main ion [amu]
constants.plasma.Zi = 1;                          % Charge number of main ion (H+) [e]
constants.plasma.Zeff = 3;                        % Typical effective charge for TCABR discharges
constants.plasma.kvPS = -1.83;                    % Pfirsch-Schlüter regime coefficient

% --- Analysis Settings ---
constants.analysis.num_radial_points = 1000;      % Number of points in the fine radial grid
constants.analysis.num_zeros_bessel = 70;         % Number of Bessel zeros for the series fit
constants.analysis.num_iterations = 5000;         % Number of iterations for bootstrap/Monte Carlo methods

% --- Default Model-Specific Parameters ---
% Store different sets of parameters here for easy switching and reproducibility.

% Set 1: Parameters from the original thesis analysis
constants.defaults.neutrals_thesis.n_H0_centre = 2e15;      % [m^-3]
constants.defaults.neutrals_thesis.n_H0_max_amp = 5e16;       % [m^-3]
constants.defaults.neutrals_thesis.r_max_H_alpha = 0.16;      % [m]
constants.defaults.neutrals_thesis.width = 0.04;          % [m]

% Set 2: Parameters derived from fitting to scaling laws
constants.defaults.neutrals_fit.n_H0_centre = 1e15;       % [m^-3]
constants.defaults.neutrals_fit.n_H0_max_amp = 2.2e16;    % [m^-3]
constants.defaults.neutrals_fit.r_max_H_alpha = 0.1564;   % [m]
constants.defaults.neutrals_fit.width = 0.022;            % [m]

% --- Active Model Parameters ---
% Select which set of defaults to use for the current run.
% To change the model, simply change the source struct here.
constants.models.neutrals = constants.defaults.neutrals_thesis; % for the validation run
% constants.models.neutrals = constants.defaults.neutrals_fit; % for neutrals Voigt fit

% --- File Paths ---
[base_dir, ~, ~] = fileparts(mfilename('fullpath'));
project_root = base_dir; % Assumes this file is in 'src' or root depending on layout. 
% Adjust if setup_constants is in utils/ or src/ and needs to point to ../data
% If setup_constants is in root:
constants.data.filename = 'experimental_profiles.txt';
constants.data.path = fullfile(project_root, 'data', constants.data.filename);

end