function [derived_profiles] = compute_derived_profiles(exp_data, ...
    temperature_results, magnetic_field, constants, r_fine)
%COMPUTE_DERIVED_PROFILES Calculates key derived physics profiles for transport analysis.
%   This function serves as a central hub for calculating essential plasma
%   profiles derived from the basic experimental and fitted data. It computes
%   parameters such as the safety factor (q), magnetic shear (s), thermal
%   velocities, various collision metrics, and standard collisionality (nu*).
%
%   Syntax:
%       derived_profiles = compute_derived_profiles(exp_data, temperature_results, magnetic_field, constants)
%
%   Inputs:
%       exp_data           - Structure with raw experimental profiles (ne, ni).
%       temperature_results - Structure with fitted temperature profiles.
%       magnetic_field     - Structure with calculated magnetic field profiles.
%       constants          - Structure with all project constants.
%
%   Output:
%       derived_profiles - A structure containing all calculated derived profiles.

fprintf('Calculating derived physics profiles...\n');

% --- 1. Extract Core Variables and Prepare Grids ---

B_tor   = magnetic_field.toroidal;
B_pol   = magnetic_field.poloidal;
a       = constants.machine.a;
R0      = constants.machine.R0;
me      = constants.physics.me;
mp      = constants.physics.mp; % Main ion is Hydrogen (proton)
e       = constants.physics.e_charge;

% Use the raw density profiles directly, as they are already on the fine grid.
% The unnecessary and erroneous interpolation step has been removed.
ne_fine = exp_data.ne_raw;
ni_fine = exp_data.ni_raw;

% Get the mean fitted temperature profiles
Ti_mean = temperature_results.profile_avg; % [eV]
% Electron temperature profile is calculated analytically, based on the
% fitted ion temperature at the edge, as done in the original script.
Te_mean = (509 - Ti_mean(end)) .* (1 - (r_fine / a).^2).^3.3 + Ti_mean(end); % [eV]


% --- 2. Safety Factor (q) and Magnetic Shear (s) ---
% Using the more accurate formula with toroidal correction from the original script.
fprintf('... calculating safety factor and magnetic shear\n');

B_pol_safe = B_pol + 1e-9; % Add small offset to avoid division by zero at the centre
% 'fq' is an empirical factor (fix factor by Severo & Novaes, 08/07/2024)
% used to refine the q-profile calculation.
fq = exp(0.21303 .* r_fine ./ a); 
q_profile = (r_fine ./ R0) .* (B_tor ./ B_pol_safe) ./ sqrt(1 - (r_fine ./ R0).^2) .* fq;

% The formula results in NaN at r=0 (0/0). We handle this by extrapolating
% from the second point, which is a common and physically reasonable approach.
q_profile(1) = q_profile(2); 

s_profile = (r_fine ./ q_profile) .* gradient(q_profile, r_fine);
s_profile(1) = 0; % Shear is zero at the magnetic axis by definition.


% --- 3. Thermal Velocities ---
% The thermal velocity is calculated from the fundamental definition: v_th = sqrt(2*k*T/m).
% This was cross-validated in the original script against the textbook formula from
% WESSON (2.4 Larmor orbits, p. 42), which is v_T ~ 4.3e5 * sqrt(T[keV]/A_i),
% confirming their equivalence after accounting for constants and unit conversions.
fprintf('... calculating thermal velocities\n');
v_th_i = sqrt(2 * e * Ti_mean / mp); % Ion thermal velocity [m/s]
v_th_e = sqrt(2 * e * Te_mean / me); % Electron thermal velocity [m/s]


% --- 4. Coulomb Logarithms and Collision Frequencies ---
% These calculations follow the standard Spitzer-Braginskii model.
fprintf('... calculating collision frequencies\n');

% Coulomb logarithms (dimensionless), formulae from WESSON (14.5 Coulomb logarithm, p. 727).
ln_lambda_ee = 14.9 - 0.5 * log(ne_fine./1e20) + log(Te_mean./1000);
ln_lambda_ei = 15.2 - 0.5 * log(ne_fine ./ 1e20) + log(Te_mean ./ 1000);
ln_lambda_ii = 17.3 - 0.5 * log(ni_fine ./ 1e20) + 1.5 * log(Ti_mean ./ 1000);

% Ion-ion collision time [s], from WESSON (14.6 Collision times, Eq. 14.6.2, p. 730).
tau_i_calc = 6.6e17 * sqrt(constants.plasma.m_ion_amu) .* (Ti_mean./1000).^(3/2) ./ (ni_fine .* ln_lambda_ii);
tau_i_15   = 6.6e17 * sqrt(constants.plasma.m_ion_amu) .* (Ti_mean./1000).^(3/2) ./ (ni_fine .* 15); % <-- ADICIONADO
tau_i_17   = 6.6e17 * sqrt(constants.plasma.m_ion_amu) .* (Ti_mean./1000).^(3/2) ./ (ni_fine .* 17); % <-- ADICIONADO

% Electron-ion collision time [s], from WESSON (14.6 Collision times, Eq. 14.6.1, p. 729).
tau_e_calc = 1.09e16 .* (Te_mean./1000).^(3/2) ./ (ni_fine .* ln_lambda_ei);
tau_e_15   = 1.09e16 .* (Te_mean./1000).^(3/2) ./ (ni_fine .* 15);
tau_e_17   = 1.09e16 .* (Te_mean./1000).^(3/2) ./ (ni_fine .* 17);


nu_i_calc = 1 ./ tau_i_calc; % Ion-ion collision frequency [s^-1]
nu_ei_calc = 1 ./ tau_e_calc; % Electron-ion collision frequency [s^-1]

nu_i_15 = 1 ./ tau_i_15; % Ion-ion collision frequency [s^-1]
nu_ei_15 = 1 ./ tau_e_15; % Electron-ion collision frequency [s^-1]

% --- 5. Collisionality (nu_*) ---
% This is the standard definition of collisionality, as defined in WESSON (14.12 Bootstrap current, p.739).
% NOTE: To maintain consistency with the original monolithic code and the results
% presented in the thesis, this calculation uses the global inverse aspect ratio
% (epsilon = a/R0) rather than the local value (epsilon_r = r/R0). This avoids
% the numerical divergence at the magnetic axis.
fprintf('... calculating collisionality\n');

R_coord = R0 * (1 + r_fine / R0);
epsilon = constants.machine.epsilon_aspect_ratio; % Use the global, constant value

% Bounce frequencies [rad/s]
omega_bounce_i = sqrt(epsilon) .* v_th_i ./ (q_profile .* R_coord);
omega_bounce_e = sqrt(epsilon) .* v_th_e ./ (q_profile .* R_coord);

% Collisionality (dimensionless)
nu_star_i = nu_i_15 ./ (epsilon .* omega_bounce_i);
nu_star_e = nu_ei_15 ./ (epsilon .* omega_bounce_e);


% --- 6. Normalised Density Gradient (R/Ln) ---
L_ne_inv = -gradient(log(ne_fine), r_fine); % Inverse density gradient scale length, 1/L_n [m^-1]
R_over_Ln = R0 .* L_ne_inv; % Normalised density gradient, R/L_n


% --- 7. Package Results into Output Structure ---
derived_profiles.r_fine     = r_fine;
derived_profiles.q_profile  = q_profile;
derived_profiles.s_profile  = s_profile;
derived_profiles.Te_profile = Te_mean;
derived_profiles.v_th_i     = v_th_i;
derived_profiles.v_th_e     = v_th_e;
derived_profiles.ln_lambda_ee = ln_lambda_ee; % <-- ADICIONADO
derived_profiles.ln_lambda_ei = ln_lambda_ei; % <-- ADICIONADO
derived_profiles.ln_lambda_ii = ln_lambda_ii; % <-- ADICIONADO
derived_profiles.tau_e_calc = tau_e_calc;     % <-- ADICIONADO
derived_profiles.tau_e_15 = tau_e_15;         % <-- ADICIONADO
derived_profiles.tau_e_17 = tau_e_17;         % <-- ADICIONADO
derived_profiles.tau_i_calc = tau_i_calc;     % <-- ADICIONADO
derived_profiles.tau_i_15 = tau_i_15;         % <-- ADICIONADO
derived_profiles.tau_i_17 = tau_i_17;         % <-- ADICIONADO
derived_profiles.nu_i_calc  = nu_i_calc;
derived_profiles.nu_ei_calc = nu_ei_calc;
derived_profiles.nu_i_15    = nu_i_15;
derived_profiles.nu_ei_15   = nu_ei_15;
derived_profiles.nu_star_i  = nu_star_i;
derived_profiles.nu_star_e  = nu_star_e;
derived_profiles.R_over_Ln  = R_over_Ln;

fprintf('Derived profiles calculated successfully.\n\n');

end