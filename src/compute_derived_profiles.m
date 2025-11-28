function [derived_profiles] = compute_derived_profiles(exp_data, ...
    temperature_results, magnetic_field, constants, r_fine)
%COMPUTE_DERIVED_PROFILES Calculates key derived physics profiles for transport analysis.
%   This function serves as a central hub for calculating essential plasma
%   profiles derived from the basic experimental and fitted data.
%
%   CORE TASKS:
%   1.  Calculates local profiles: Safety factor (q), Shear (s), Thermal velocities.
%   2.  Computes Collision Frequencies based on Spitzer-Braginskii formalism.
%   3.  Determines Local Collisionality profiles (Wesson/Bounce and Solomon/Drift).
%   4.  **NEW:** Calculates Volume-Averaged quantities (<Te>, <ne>) and 
%       Global Representative Collisionality to serve as robust inputs for 
%       multi-machine scaling laws (avoiding edge singularities).
%
%   Syntax:
%       derived_profiles = compute_derived_profiles(exp_data, temperature_results, ...
%                                                   magnetic_field, constants, r_fine)
%
%   Inputs:
%       exp_data           - Structure with raw experimental profiles (ne, ni).
%       temperature_results - Structure with fitted temperature profiles.
%       magnetic_field     - Structure with calculated magnetic field profiles.
%       constants          - Structure with all project constants.
%       r_fine             - The high-resolution radial grid [m].
%
%   Output:
%       derived_profiles - A structure containing all calculated derived profiles
%                          and global parameters.

    arguments
        exp_data (1,1) struct
        temperature_results (1,1) struct
        magnetic_field (1,1) struct
        constants (1,1) struct
        r_fine (:,1) {mustBeNumeric, mustBeReal, mustBeFinite}
    end

    fprintf('Calculating derived physics profiles...\n');

    % =========================================================================
    % 1. EXTRACT CORE VARIABLES AND PREPARE GRIDS
    % =========================================================================
    B_tor   = magnetic_field.toroidal;
    B_pol   = magnetic_field.poloidal;
    a       = constants.machine.a;
    R0      = constants.machine.R0;
    me      = constants.physics.me;
    mp      = constants.physics.mp; % Main ion is Hydrogen (proton)
    e       = constants.physics.e_charge;
    
    % Experimental Density Profiles (already interpolated to r_fine if loaded correctly,
    % but ensuring consistency here).
    ne_fine = exp_data.ne_raw;
    ni_fine = exp_data.ni_raw;
    
    % Temperature Profiles
    % Ti: From canonical fit (bootstrap mean)
    Ti_mean = temperature_results.profile_avg; % [eV]
    
    % Te: Calculated analytically based on Ti fit shape (Assumption from thesis methodology)
    % T_e(r) profile shape follows T_i but scaled.
    Te_mean = (509 - Ti_mean(end)) .* (1 - (r_fine / a).^2).^3.3 + Ti_mean(end); % [eV]

    % =========================================================================
    % 2. SAFETY FACTOR (q) AND MAGNETIC SHEAR (s)
    % =========================================================================
    fprintf('... calculating safety factor and magnetic shear\n');
    q_and_shear = struct();
    
    B_pol_safe = B_pol + 1e-9; % Avoid division by zero at axis
    
    % 'fq' is an empirical geometric factor for TCABR (Severo & Novaes, 2024)
    fq = exp(0.21303 .* r_fine ./ a); 
    
    % Cylindrical approximation with toroidal correction
    q_profile = (r_fine ./ R0) .* (B_tor ./ B_pol_safe) ./ sqrt(1 - (r_fine ./ R0).^2) .* fq;
    
    % Handle singularity at r=0 by extrapolating from the first valid point
    q_profile(1) = q_profile(2); 

    % Magnetic Shear: s = (r/q) * (dq/dr)
    s_profile_local = (r_fine ./ q_profile) .* gradient(q_profile, r_fine);
    s_profile_local(1) = 0; % Enforce zero shear at axis

    q_and_shear.q_profile = q_profile;
    q_and_shear.s_profile = s_profile_local;

    % =========================================================================
    % 3. THERMAL VELOCITIES
    % =========================================================================
    fprintf('... calculating thermal velocities\n');
    thermal_velocities = struct();
    thermal_velocities.v_th_i = sqrt(2 * e * Ti_mean / mp); % [m/s]
    thermal_velocities.v_th_e = sqrt(2 * e * Te_mean / me); % [m/s]

    % =========================================================================
    % 4. COLLISIONS (FREQUENCIES AND TIMES)
    % =========================================================================
    % Based on Wesson (Tokamaks, 3rd Ed.), Chapter 14.
    fprintf('... calculating collision metrics\n');
    collisions = struct();
    
    % Coulomb Logarithms
    collisions.ln_lambda_ee = 14.9 - 0.5 * log(ne_fine./1e20) + log(Te_mean./1000);
    collisions.ln_lambda_ei = 15.2 - 0.5 * log(ne_fine ./ 1e20) + log(Te_mean ./ 1000);
    collisions.ln_lambda_ii = 17.3 - 0.5 * log(ni_fine ./ 1e20) + 1.5 * log(Ti_mean ./ 1000);

    % Collision Times [s]
    % Ion-Ion
    tau_i_calc = 6.6e17 * sqrt(constants.plasma.m_ion_amu) .* (Ti_mean./1000).^(3/2) ./ (ni_fine .* collisions.ln_lambda_ii);
    % Fixed Lambda approximations for comparison/legacy plots
    tau_i_15   = 6.6e17 * sqrt(constants.plasma.m_ion_amu) .* (Ti_mean./1000).^(3/2) ./ (ni_fine .* 15); 
    tau_i_17   = 6.6e17 * sqrt(constants.plasma.m_ion_amu) .* (Ti_mean./1000).^(3/2) ./ (ni_fine .* 17); 

    % Electron-Ion
    tau_e_calc = 1.09e16 .* (Te_mean./1000).^(3/2) ./ (ni_fine .* collisions.ln_lambda_ei);
    tau_e_15   = 1.09e16 .* (Te_mean./1000).^(3/2) ./ (ni_fine .* 15);
    tau_e_17   = 1.09e16 .* (Te_mean./1000).^(3/2) ./ (ni_fine .* 17);

    % Store Times
    collisions.tau_i_calc = tau_i_calc; collisions.tau_i_15 = tau_i_15; collisions.tau_i_17 = tau_i_17;
    collisions.tau_e_calc = tau_e_calc; collisions.tau_e_15 = tau_e_15; collisions.tau_e_17 = tau_e_17;

    % Collision Frequencies [s^-1]
    collisions.nu_i_calc  = 1 ./ tau_i_calc;
    collisions.nu_ei_calc = 1 ./ tau_e_calc;
    collisions.nu_i_15    = 1 ./ tau_i_15;
    collisions.nu_ei_15   = 1 ./ tau_e_15;

    % =========================================================================
    % 5. COLLISIONALITY (LOCAL & GLOBAL)
    % =========================================================================
    fprintf('... calculating collisionality (Local & Global)\n');
    collisionality = struct();
    
    R_coord = R0 * (1 + r_fine / R0);
    epsilon_global = constants.machine.epsilon_aspect_ratio;

    % --- 5a. Local Collisionality (Wesson / Bounce) ---
    % nu* = nu_eff / omega_bounce
    omega_bounce_i = sqrt(epsilon_global) .* thermal_velocities.v_th_i ./ (q_and_shear.q_profile .* R_coord);
    omega_bounce_e = sqrt(epsilon_global) .* thermal_velocities.v_th_e ./ (q_and_shear.q_profile .* R_coord);

    collisionality.nu_star_i_wesson = collisions.nu_i_15 ./ (epsilon_global .* omega_bounce_i);
    collisionality.nu_star_e_wesson = collisions.nu_ei_15 ./ (epsilon_global .* omega_bounce_e);
    
    % Aliases
    collisionality.nu_star_i = collisionality.nu_star_i_wesson;
    collisionality.nu_star_e = collisionality.nu_star_e_wesson;

    % --- 5b. Local Collisionality (Solomon / Maslov Drift) ---
    % nu* = nu_ei / omega_De (Curvature Drift Frequency)
    % Used for local profile diagnostics, but known to diverge at cold edge.
    c_s = sqrt(e * Te_mean ./ mp); 
    k_perp_rho_s = sqrt(0.1); 
    omega_De = 2 * k_perp_rho_s .* c_s ./ R_coord;
    
    collisionality.nu_star_e_solomon_local = collisions.nu_ei_15 ./ omega_De;

    % --- 5c. GLOBAL REPRESENTATIVE COLLISIONALITY (The "Correction") ---
    % Consistent with Maslov (2009) and Solomon (2010) methodology:
    % Use Volume-Averaged parameters <Te>, <ne> to define the discharge regime.
    % This avoids edge singularities and provides a robust scalar input for scaling laws.
    
    % 1. Calculate Volume Averages <Y> = int(Y * r dr) / int(r dr)
    % Differential volume element dV ~ r (cylindrical approx is sufficient for weighting)
    r_weight = r_fine; 
    
    % Helper function for trapezoidal integration
    integrate_vol = @(y) trapz(r_fine, y .* r_weight);
    vol_norm = integrate_vol(ones(size(r_fine)));
    
    Te_vol_avg = integrate_vol(Te_mean) / vol_norm;
    ne_vol_avg = integrate_vol(ne_fine) / vol_norm;
    
    fprintf('      <Te> (Vol. Avg): %.2f eV\n', Te_vol_avg);
    fprintf('      <ne> (Vol. Avg): %.2e m^-3\n', ne_vol_avg);
    
    % 2. Calculate Global Collision Frequency (at <Te>, <ne>)
    % nu_ei_global ~ ne * Te^-1.5
    ln_lambda_global = 15.2 - 0.5 * log(ne_vol_avg/1e20) + log(Te_vol_avg/1000);
    tau_e_global = 1.09e16 * (Te_vol_avg/1000)^(3/2) / (ne_vol_avg * ln_lambda_global);
    nu_ei_global = 1 / tau_e_global;
    
    % 3. Calculate Global Drift Frequency (at <Te>, R0)
    % omega_d ~ Te^0.5 / R0
    c_s_global = sqrt(e * Te_vol_avg / mp);
    omega_De_global = 2 * sqrt(0.1) * c_s_global / R0;
    
    % 4. Global Solomon Collisionality
    nu_star_e_global_solomon = nu_ei_global / omega_De_global;
    
    fprintf('      Global Solomon Collisionality (nu*_global): %.4f\n', nu_star_e_global_solomon);
    
    % Save to structure
    collisionality.global.Te_avg = Te_vol_avg;
    collisionality.global.ne_avg = ne_vol_avg;
    collisionality.global.nu_star_e_solomon = nu_star_e_global_solomon;
    
    % For compatibility, we create a "profile" that is constant, 
    % but the advanced logic will prefer the scalar value where appropriate.
    collisionality.nu_star_e_solomon = collisionality.nu_star_e_solomon_local; % Default to local for plots

    % =========================================================================
    % 6. GRADIENTS
    % =========================================================================
    gradients = struct();
    L_ne_inv = -gradient(log(ne_fine), r_fine); 
    gradients.R_over_Ln = R0 .* L_ne_inv;

    % =========================================================================
    % 7. PACKAGE RESULTS
    % =========================================================================
    derived_profiles.r_fine     = r_fine;
    derived_profiles.q_profile  = q_profile;
    derived_profiles.s_profile  = s_profile_local;
    derived_profiles.Te_profile = Te_mean;
    derived_profiles.v_th_i     = thermal_velocities.v_th_i;
    derived_profiles.v_th_e     = thermal_velocities.v_th_e;
    
    % Flatten collision fields
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

    % Structure fields
    derived_profiles.thermal_velocities = thermal_velocities;
    derived_profiles.collisions = collisions;
    derived_profiles.collisionality = collisionality;
    
    % Shortcut fields
    derived_profiles.nu_star_i = collisionality.nu_star_i;
    derived_profiles.nu_star_e = collisionality.nu_star_e;
    derived_profiles.nu_star_e_solomon = collisionality.nu_star_e_solomon; % Local profile
    
    % Global Params Access
    derived_profiles.global_params = collisionality.global;

    derived_profiles.R_over_Ln  = gradients.R_over_Ln;

    fprintf('Derived profiles calculated successfully.\n\n');
end