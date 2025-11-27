function [derived_profiles] = compute_derived_profiles(exp_data, ...
    temperature_results, magnetic_field, constants, r_fine)
%COMPUTE_DERIVED_PROFILES Calculates key derived physics profiles for transport analysis.
%   This function serves as a central hub for calculating essential plasma
%   profiles derived from the basic experimental and fitted data. It computes
%   parameters such as the safety factor (q), magnetic shear (s), thermal
%   velocities, various collision metrics, and standard collisionality (nu*).
%
%   Syntax:
%       derived_profiles = compute_derived_profiles(exp_data, temperature_results, magnetic_field, constants, r_fine)
%
%   Inputs:
%       exp_data           - Structure with raw experimental profiles (ne, ni).
%       temperature_results - Structure with fitted temperature profiles.
%       magnetic_field     - Structure with calculated magnetic field profiles.
%       constants          - Structure with all project constants.
%       r_fine             - Radial grid.
%
%   Output:
%       derived_profiles - A structure containing all calculated derived profiles.

    arguments
        exp_data (1,1) struct
        temperature_results (1,1) struct
        magnetic_field (1,1) struct
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
    
    % Get the mean fitted temperature profiles
    Ti_mean = temperature_results.profile_avg; % [eV]
    
    % Electron temperature profile is calculated analytically, based on the
    % fitted ion temperature at the edge, as done in the original script.
    Te_mean = (509 - Ti_mean(end)) .* (1 - (r_fine / a).^2).^3.3 + Ti_mean(end); % [eV]

    % --- 2. Safety Factor (q) and Magnetic Shear (s) ---
    fprintf('... calculating safety factor and magnetic shear\n');
    q_and_shear = struct();
    
    B_pol_safe = B_pol + 1e-9; % Add small offset to avoid division by zero at the centre
    % 'fq' is an empirical factor (fix factor by Severo & Novaes, 08/07/2024)
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

    % --- 4. Coulomb Logarithms and Collision Frequencies ---
    fprintf('... calculating collision frequencies\n');
    collisions = struct();
    
    collisions.ln_lambda_ee = 14.9 - 0.5 * log(ne_fine./1e20) + log(Te_mean./1000);
    collisions.ln_lambda_ei = 15.2 - 0.5 * log(ne_fine ./ 1e20) + log(Te_mean ./ 1000);
    collisions.ln_lambda_ii = 17.3 - 0.5 * log(ni_fine ./ 1e20) + 1.5 * log(Ti_mean ./ 1000);

    % Ion-ion collision time [s]
    tau_i_calc = 6.6e17 * sqrt(constants.plasma.m_ion_amu) .* (Ti_mean./1000).^(3/2) ./ (ni_fine .* collisions.ln_lambda_ii);
    tau_i_15   = 6.6e17 * sqrt(constants.plasma.m_ion_amu) .* (Ti_mean./1000).^(3/2) ./ (ni_fine .* 15); 
    tau_i_17   = 6.6e17 * sqrt(constants.plasma.m_ion_amu) .* (Ti_mean./1000).^(3/2) ./ (ni_fine .* 17); 

    % Electron-ion collision time [s]
    tau_e_calc = 1.09e16 .* (Te_mean./1000).^(3/2) ./ (ni_fine .* collisions.ln_lambda_ei);
    tau_e_15   = 1.09e16 .* (Te_mean./1000).^(3/2) ./ (ni_fine .* 15);
    tau_e_17   = 1.09e16 .* (Te_mean./1000).^(3/2) ./ (ni_fine .* 17);

    % Save tau variables to the struct so they can be accessed later
    collisions.tau_i_calc = tau_i_calc;
    collisions.tau_i_15   = tau_i_15;
    collisions.tau_i_17   = tau_i_17;
    collisions.tau_e_calc = tau_e_calc;
    collisions.tau_e_15   = tau_e_15;
    collisions.tau_e_17   = tau_e_17;

    collisions.nu_i_calc  = 1 ./ tau_i_calc;
    collisions.nu_ei_calc = 1 ./ tau_e_calc;
    collisions.nu_i_15    = 1 ./ tau_i_15;
    collisions.nu_ei_15   = 1 ./ tau_e_15;

    % --- 5. Collisionality (nu_*) ---
    fprintf('... calculating collisionality\n');
    collisionality = struct();
    
    R_coord = R0 * (1 + r_fine / R0);
    epsilon = constants.machine.epsilon_aspect_ratio; % Global epsilon

    % 5a. Wesson Collisionality (Standard/Bounce)
    omega_bounce_i = sqrt(epsilon) .* thermal_velocities.v_th_i ./ (q_and_shear.q_profile .* R_coord);
    omega_bounce_e = sqrt(epsilon) .* thermal_velocities.v_th_e ./ (q_and_shear.q_profile .* R_coord);

    collisionality.nu_star_i_wesson = collisions.nu_i_15 ./ (epsilon .* omega_bounce_i);
    collisionality.nu_star_e_wesson = collisions.nu_ei_15 ./ (epsilon .* omega_bounce_e);
    
    % Aliases for backward compatibility
    collisionality.nu_star_i = collisionality.nu_star_i_wesson;
    collisionality.nu_star_e = collisionality.nu_star_e_wesson;

    % 5b. Solomon/Maslov Collisionality (Drift Frequency)
    % Definition: nu*_e = nu_ei / omega_De
    % omega_De = 2 * k_perp * rho_s * c_s / R
    % with k_perp*rho_s ~ sqrt(0.1) [Maslov 2009]
    c_s = sqrt(e * Te_mean ./ mp); % Ion sound speed (approx)
    
    k_perp_rho_s = sqrt(0.1); 
    omega_De = 2 * k_perp_rho_s .* c_s ./ R_coord;
    
    collisionality.nu_star_e_solomon = collisions.nu_ei_15 ./ omega_De;

    % --- 6. Normalised Density Gradient (R/Ln) ---
    gradients = struct();
    L_ne_inv = -gradient(log(ne_fine), r_fine); 
    gradients.R_over_Ln = R0 .* L_ne_inv;

    % --- 7. Package Results into Output Structure ---
    derived_profiles.r_fine     = r_fine;
    derived_profiles.q_profile  = q_and_shear.q_profile;
    derived_profiles.s_profile  = q_and_shear.s_profile;
    derived_profiles.Te_profile = Te_mean;
    derived_profiles.v_th_i     = thermal_velocities.v_th_i;
    derived_profiles.v_th_e     = thermal_velocities.v_th_e;
    derived_profiles.collisions = collisions;
    derived_profiles.thermal_velocities = thermal_velocities;
    
    % Flatten collision fields for compatibility
    derived_profiles.ln_lambda_ee = collisions.ln_lambda_ee;
    derived_profiles.ln_lambda_ei = collisions.ln_lambda_ei;
    derived_profiles.ln_lambda_ii = collisions.ln_lambda_ii;
    derived_profiles.tau_e_calc = collisions.tau_e_calc;
    derived_profiles.tau_e_15 = collisions.tau_e_15;
    derived_profiles.tau_e_17 = collisions.tau_e_17;
    derived_profiles.tau_i_calc = collisions.tau_i_calc;
    derived_profiles.tau_i_15 = collisions.tau_i_15;
    derived_profiles.tau_i_17 = collisions.tau_i_17;
    derived_profiles.nu_i_calc = collisions.nu_i_calc;
    derived_profiles.nu_ei_calc = collisions.nu_ei_calc;
    derived_profiles.nu_i_15 = collisions.nu_i_15;
    derived_profiles.nu_ei_15 = collisions.nu_ei_15;
    
    derived_profiles.collisionality = collisionality;
    derived_profiles.nu_star_i = collisionality.nu_star_i;
    derived_profiles.nu_star_e = collisionality.nu_star_e;
    derived_profiles.nu_star_e_solomon = collisionality.nu_star_e_solomon;
    derived_profiles.nu_star_e_wesson = collisionality.nu_star_e_wesson;

    derived_profiles.R_over_Ln  = gradients.R_over_Ln;

    fprintf('Derived profiles calculated successfully.\n\n');
end