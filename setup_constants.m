function constants = setup_constants()
%SETUP_CONSTANTS Defines and returns a structure with all project constants.
%   This function organises all physical constants, machine parameters, and
%   analysis settings into a nested structure for clarity and easy access.
%
%   Output:
%       constants - A structure containing all defined constants.

% --- Machine and Discharge Parameters ---
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
constants.physics.kB = 1.380649E-23;              % Boltzmann constant [J/K]
constants.physics.c_light = 299792458;           % Speed of light [m/s]

% --- Plasma Species and Analysis Parameters ---
% NOTE: Throughout the code, 'ion' refers to the main hydrogen species (H+).
% The toroidal velocity is measured from C6+ impurities, assuming V_phi_C6+ ~ V_phi_H+.
constants.plasma.m_ion = constants.physics.mp;      % Mass of main ion (Hydrogen) [kg]
constants.plasma.m_ion_amu = 1.00784;             % Atomic mass of main ion [amu]
constants.plasma.m_impurity = 1.9926465E-26;      % Mass of impurity (Carbon-12) [kg]
constants.plasma.Zi = 1;                          % Charge number of main ion (H+) [e]
% constants.plasma.Z_impurity = 6;                  % Charge number of impurity (C6+) [e]
constants.plasma.Zeff = 3;                        % Typical effective charge
constants.plasma.coulomb_log = 17;                % Typical Coulomb logarithm
constants.plasma.kvPS = -1.83;                    % Pfirsch-Schlüter regime coefficient

% --- Analysis Settings ---
constants.analysis.num_radial_points = 1000;      % Number of points in the fine radial grid
constants.analysis.num_zeros_bessel = 70;         % Number of Bessel zeros for the series fit
constants.analysis.num_iterations = 5000;         % Number of iterations for bootstrap/Monte Carlo methods
constants.analysis.generate_final_plots = true;

% --- Model-Specific Parameters ---
% Parameters for the Gaussian neutral density profile model (n_H0)
% Based on thesis section 3.4, referencing Kantor et al. (2007, 2012) and Severo et al. (2021).
constants.models.neutrals.n_H0_centre = 2e15;       % Density at the centre [m^-3] [cite: 320]
constants.models.neutrals.n_H0_max_amp = 5e16;      % Maximum amplitude of the Gaussian peak [m^-3] [cite: 352]
constants.models.neutrals.r_max_H_alpha = 0.16;     % Radial position of max H-alpha emission [m]
constants.models.neutrals.width = 0.04;           % Gaussian width (penetration depth) [m] [cite: 353]

% Based on fitting to scaling laws
% constants.models.neutrals.n_H0_centre = 1e15;       % Density at the centre [m^-3] [cite: 320]
% constants.models.neutrals.n_H0_max_amp = 2.2e16;      % Maximum amplitude of the Gaussian peak [m^-3] [cite: 352]
% constants.models.neutrals.r_max_H_alpha = 0.1564;     % Radial position of max H-alpha emission [m]
% constants.models.neutrals.width = 0.022;           % Gaussian width (penetration depth) [m] [cite: 353]

% --- File Paths ---
constants.data.filename = 'experimental_profiles.txt';
constants.data.path = fullfile('data', constants.data.filename);

end