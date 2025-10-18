function constants = setup_constants()
%% SETUP_CONSTANTS
% Defines and returns a structure containing all physical and geometrical
% tokamak constants used in the analysis, as specified by the new
% MAIN_TOKAMAK_ANALYSIS.m script.
%
% Returns:
%   constants - A structure containing all defined constants.

% Geometrical Parameters of the Tokamak
constants.a = 0.18;                       % Minor radius [m]
constants.Rcentr = 0.61;                  % Major radius [m]
constants.epsilon_aspect_ratio = constants.a / constants.Rcentr; % Inverse aspect ratio

% Fundamental Physical Constants (SI units)
constants.e_charge = 1.60217663E-19;      % Electron electric charge [C]
constants.me_electron = 9.1093837E-31;    % Electron mass [kg]
constants.mp_proton = 1.672622E-27;       % Mass of proton [kg]
constants.m_ion = 1.9926465E-26;          % Mass of Impurity/main ion [kg] (e.g., carbon, or Deuterium if it was 3.3435860E-27)
constants.m_ion_amu = 12.011;             % Mass of main ion in amu (e.g., carbon: 12.011, for Deuterium: 2.014)
constants.c_light = 299792458;            % Speed of light in vacuum [m/s]
constants.mu0_SI = 4 * pi * 1E-7;         % Permeability of free space [N/A^2 or H/m]
constants.kB_Boltzmann = 1.380649E-23;    % Boltzmann constant [J/K]
constants.epsilon0_SI = 8.854187812E-12;  % Vacuum permittivity [C^2*s^2*m^-3*kg^-1 or F/m]

% Plasma-Specific Constants
constants.Lamb = 17;                      % Coulomb logarithm (ln\Lambda ~ 17), a typical value.
constants.Zi = 1;                         % Ion charge number [e] (e.g., 1 for Hydrogen/Deuterium).
constants.ZI = 6;                         % Impurity electric charge (for C6+)
constants.Zeff = 3;                       % Effective charge of plasma.
constants.Ip = 90E3;                      % Plasma current [A] (as per original MAIN_CODE_OK.m).
constants.Bcentr = 1.07;                  % Main axis magnetic field [T] (as per original MAIN_CODE_OK.m).

% Specific coefficients from original code (if not derived elsewhere)
constants.kvPS = -1.83;                   % Collisionality Pfirsch-Schlüter regime coefficient.

% Ensure units are consistent for later use (e.g., converting eV to J)
constants.eV_to_Joule = constants.e_charge;
constants.Joule_to_eV = 1 / constants.e_charge;

% --- File Paths ---
constants.data.filename = 'experimental_profiles.txt';
constants.data.path = fullfile('data', constants.data.filename);

end