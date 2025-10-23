function [derived_profiles] = compute_derived_profiles(exp_data, ...
    temperature_results, magnetic_field, constants, r_fine)
%COMPUTE_DERIVED_PROFILES Calculates key derived physics profiles for transport analysis.
%   This function is a central hub for calculating essential plasma profiles
%   derived from experimental and fitted data. It computes quantities such as:
%   - Safety factor (q) and magnetic shear (s)
%   - Thermal velocities (ion and electron)
%   - Coulomb logarithms and collision frequencies/times
%   - Standard collisionality (nu*)
%   - Normalised density gradient (R/Ln)
%
%   Syntax:
%       derived_profiles = compute_derived_profiles(exp_data, ...
%           temperature_results, magnetic_field, constants, r_fine)
%
%   Output:
%       derived_profiles - A structure with nested fields containing all profiles.

    arguments
        exp_data (1,1) struct {mustContainFields(exp_data, {'ne_raw', 'ni_raw'})}
        temperature_results (1,1) struct {mustContainFields(temperature_results, {'profile_avg'})}
        magnetic_field (1,1) struct {mustContainFields(magnetic_field, {'toroidal', 'poloidal'})}
        constants (1,1) struct
        r_fine (:,1) {mustBeNumeric, mustBeReal, mustBeFinite}
    end

    fprintf('Calculating derived physics profiles...\n');

    % --- 1. Extract Core Variables and Prepare Grids ---
    B_tor   = magnetic_field.toroidal;
    B_pol   = magnetic_field.poloidal;
    a       = constants.machine.a;
    R0      = constants.machine.R0;
    me      = constants.physics.me;
    mp      = constants.physics.mp; % Main ion is Hydrogen (proton)
    e       = constants.physics.e_charge;

    ne_fine = exp_data.ne_raw;
    ni_fine = exp_data.ni_raw;
    Ti_mean = temperature_results.profile_avg; % [eV]

    % Electron temperature profile is calculated analytically, as per original analysis.
    Te_mean = (509 - Ti_mean(end)) .* (1 - (r_fine / a).^2).^3.3 + Ti_mean(end); % [eV]

    % --- 2. Safety Factor (q) and Magnetic Shear (s) ---
    fprintf('... calculating safety factor and magnetic shear\n');
    q_and_shear = struct();
    
    % Use small offset in B_pol to avoid division by zero at the magnetic axis
    B_pol_safe = B_pol + 1e-9;
    
    % 'fq' is an empirical factor (Severo & Novaes, 08/07/2024) used to
    % refine the q-profile calculation based on experimental constraints.
    fq = exp(0.21303 .* r_fine ./ a);
    q_profile = (r_fine ./ R0) .* (B_tor ./ B_pol_safe) ./ sqrt(1 - (r_fine ./ R0).^2) .* fq;
    
    % Handle NaN at r=0 (0/0) by extrapolating from the second point.
    q_profile(1) = q_profile(2);

    q_and_shear.q_profile = q_profile;
    q_and_shear.s_profile = (r_fine ./ q_profile) .* gradient(q_profile, r_fine);
    q_and_shear.s_profile(1) = 0; % Shear is zero at the axis by definition.

    % --- 3. Thermal Velocities ---
    fprintf('... calculating thermal velocities\n');
    thermal_velocities = struct();
    thermal_velocities.v_th_i = sqrt(2 * e * Ti_mean / mp); % [m/s]
    thermal_velocities.v_th_e = sqrt(2 * e * Te_mean / me); % [m/s]

    % --- 4. Collisions (Logarithms, Times, Frequencies) ---
    fprintf('... calculating collision metrics\n');
    collisions = struct();
    
    % Coulomb logarithms (Wesson, 14.5, p. 727)
    collisions.ln_lambda_ee = 14.9 - 0.5 * log(ne_fine./1e20) + log(Te_mean./1000);
    collisions.ln_lambda_ei = 15.2 - 0.5 * log(ne_fine ./ 1e20) + log(Te_mean ./ 1000);
    collisions.ln_lambda_ii = 17.3 - 0.5 * log(ni_fine ./ 1e20) + 1.5 * log(Ti_mean ./ 1000);

    % Electron-ion collision time [s] (Wesson, 14.6.1, p. 729)
    collisions.tau_e_calc = 1.09e16 .* (Te_mean./1000).^(3/2) ./ (ni_fine .* collisions.ln_lambda_ei);
    collisions.tau_e_15   = 1.09e16 .* (Te_mean./1000).^(3/2) ./ (ni_fine .* 15);
    collisions.tau_e_17   = 1.09e16 .* (Te_mean./1000).^(3/2) ./ (ni_fine .* 17);

    % Ion-ion collision time [s] (Wesson, 14.6.2, p. 730)
    collisions.tau_i_calc = 6.6e17 * sqrt(constants.plasma.m_ion_amu) .* (Ti_mean./1000).^(3/2) ./ (ni_fine .* collisions.ln_lambda_ii);
    collisions.tau_i_15   = 6.6e17 * sqrt(constants.plasma.m_ion_amu) .* (Ti_mean./1000).^(3/2) ./ (ni_fine .* 15);
    collisions.tau_i_17   = 6.6e17 * sqrt(constants.plasma.m_ion_amu) .* (Ti_mean./1000).^(3/2) ./ (ni_fine .* 17);

    % Frequencies are inverse of times
    collisions.nu_ei_calc = 1 ./ collisions.tau_e_calc;
    collisions.nu_ei_15   = 1 ./ collisions.tau_e_15;
    collisions.nu_i_calc  = 1 ./ collisions.tau_i_calc;
    collisions.nu_i_15    = 1 ./ collisions.tau_i_15;

    % --- 5. Collisionality (nu_*) ---
    fprintf('... calculating collisionality\n');
    collisionality = struct();
    
    % NOTE: Using the global inverse aspect ratio (epsilon = a/R0) instead of
    % the local value (r/R0) for consistency with the thesis and to avoid
    % numerical divergence at the magnetic axis.
    epsilon = constants.machine.epsilon_aspect_ratio;
    R_coord = R0 * (1 + r_fine / R0);

    omega_bounce_i = sqrt(epsilon) .* thermal_velocities.v_th_i ./ (q_and_shear.q_profile .* R_coord);
    omega_bounce_e = sqrt(epsilon) .* thermal_velocities.v_th_e ./ (q_and_shear.q_profile .* R_coord);

    collisionality.nu_star_i = collisions.nu_i_15 ./ (epsilon .* omega_bounce_i);
    collisionality.nu_star_e = collisions.nu_ei_15 ./ (epsilon .* omega_bounce_e);
    
    % --- 5b. Solomon/Maslov Collisionality ---
    % This definition is required for a consistent comparison with the Solomon et al. (2010)
    % scaling laws, which are based on the effective collisionality defined in Maslov et al. (2009).
    % It is defined as the ratio of the electron-ion collision frequency
    % to the curvature drift frequency.
    % Formula: v_eff = 10e-14 * Z_eff * R * n_e / T_e^2  (with T_e in eV).
    % The equivalent definition is v_ei / w_De, where w_De ~ T_e / (R * B).
    % This implementation uses local profile values for R, n_e, and T_e.
    % NOTE: Solomon uses Te in keV, and Maslov's formula uses <Te> in eV.
    % Maslov's formula v_eff = 10^-14 * Z_eff * R * n_e / T_e^2 (with Te in eV)
    % is the one we will implement for the radial profiles.
    fprintf('... calculating Solomon/Maslov collisionality\n');
    collisionality.nu_star_e_solomon = 1e-14 * constants.plasma.Zeff .* R_coord .* ne_fine ./ (Te_mean.^2);

    % --- 6. Normalised Density Gradient (R/Ln) ---
    gradients = struct();
    L_ne_inv = -gradient(log(ne_fine), r_fine); % 1/L_n [m^-1]
    gradients.R_over_Ln = R0 .* L_ne_inv;

    % --- 7. Package Results into Output Structure ---
    derived_profiles.r_fine = r_fine;
    derived_profiles.Te_profile = Te_mean;
    derived_profiles.q_and_shear = q_and_shear;
    derived_profiles.thermal_velocities = thermal_velocities;
    derived_profiles.collisions = collisions;
    derived_profiles.collisionality = collisionality;
    derived_profiles.gradients = gradients;

    fprintf('Derived profiles calculated successfully.\n\n');
end