function [derived_profiles] = compute_derived_profiles(exp_data, ...
    temperature_results, magnetic_field, constants, r_fine)
%COMPUTE_DERIVED_PROFILES Calculates key derived physics profiles for transport analysis.
%   ... (docstring) ...

    arguments
        exp_data (1,1) struct
        temperature_results (1,1) struct
        magnetic_field (1,1) struct
        constants (1,1) struct
        r_fine (:,1) {mustBeNumeric, mustBeReal, mustBeFinite}
    end

    fprintf('Calculating derived physics profiles...\n');

    % --- 1. Input Validation ---
    if ~isfield(exp_data, 'ne_raw') || ~isfield(exp_data, 'ni_raw')
        error('Input struct ''exp_data'' is missing required fields: ''ne_raw'', ''ni_raw''.');
    end
    if ~isfield(temperature_results, 'profile_avg')
        error('Input struct ''temperature_results'' is missing required field: ''profile_avg''.');
    end
    if ~isfield(magnetic_field, 'toroidal') || ~isfield(magnetic_field, 'poloidal')
        error('Input struct ''magnetic_field'' is missing required fields: ''toroidal'', ''poloidal''.');
    end

    % --- 2. Extract Core Variables and Prepare Grids ---
    B_tor   = magnetic_field.toroidal;
    B_pol   = magnetic_field.poloidal;
    a       = constants.machine.a;
    R0      = constants.machine.R0;
    me      = constants.physics.me;
    mp      = constants.physics.mp;
    e       = constants.physics.e_charge;
    
    ne_fine = exp_data.ne_raw;
    ni_fine = exp_data.ni_raw;
    Ti_mean = temperature_results.profile_avg;
    
    % Electron temperature profile is calculated analytically
    Te_mean = (509 - Ti_mean(end)) .* (1 - (r_fine / a).^2).^3.3 + Ti_mean(end);

    % --- 3. Safety Factor (q) and Magnetic Shear (s) ---
    fprintf('... calculating safety factor and magnetic shear\n');
    q_and_shear = struct();
    
    B_pol_safe = B_pol + 1e-9;
    fq = exp(0.21303 .* r_fine ./ a); 
    q_profile = (r_fine ./ R0) .* (B_tor ./ B_pol_safe) ./ sqrt(1 - (r_fine ./ R0).^2) .* fq;
    q_profile(1) = q_profile(2); 

    q_and_shear.q_profile = q_profile;
    q_and_shear.s_profile = (r_fine ./ q_profile) .* gradient(q_profile, r_fine);
    q_and_shear.s_profile(1) = 0;

    % --- 4. Thermal Velocities ---
    fprintf('... calculating thermal velocities\n');
    thermal_velocities = struct();
    thermal_velocities.v_th_i = sqrt(2 * e * Ti_mean / mp);
    thermal_velocities.v_th_e = sqrt(2 * e * Te_mean / me);

    % --- 5. Coulomb Logarithms and Collision Frequencies ---
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

    % **CORRECTION: Save tau variables to the struct so they can be accessed later**
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

    % --- 6. Collisionality (nu_*) ---
    fprintf('... calculating collisionality\n');
    collisionality = struct();
    
    R_coord = R0 * (1 + r_fine / R0);
    epsilon = constants.machine.epsilon_aspect_ratio; 

    omega_bounce_i = sqrt(epsilon) .* thermal_velocities.v_th_i ./ (q_and_shear.q_profile .* R_coord);
    omega_bounce_e = sqrt(epsilon) .* thermal_velocities.v_th_e ./ (q_and_shear.q_profile .* R_coord);

    % Wesson collisionality
    collisionality.nu_star_i_wesson = collisions.nu_i_15 ./ (epsilon .* omega_bounce_i);
    collisionality.nu_star_e_wesson = collisions.nu_ei_15 ./ (epsilon .* omega_bounce_e);
    
    % Generic aliases
    collisionality.nu_star_i = collisionality.nu_star_i_wesson;
    collisionality.nu_star_e = collisionality.nu_star_e_wesson;

    % --- 7. Normalised Density Gradient (R/Ln) ---
    gradients = struct();
    L_ne_inv = -gradient(log(ne_fine), r_fine); 
    gradients.R_over_Ln = R0 .* L_ne_inv;

    % --- 8. Package Results into Output Structure ---
    derived_profiles.r_fine      = r_fine;
    derived_profiles.q_and_shear = q_and_shear;
    derived_profiles.thermal_velocities = thermal_velocities;
    derived_profiles.collisions  = collisions;
    derived_profiles.collisionality = collisionality;
    derived_profiles.gradients   = gradients;
    
    % Flattening structure for plotting/compatibility
    derived_profiles.q_profile = q_and_shear.q_profile;
    derived_profiles.s_profile = q_and_shear.s_profile;
    derived_profiles.Te_profile = Te_mean;
    derived_profiles.v_th_i = thermal_velocities.v_th_i;
    derived_profiles.v_th_e = thermal_velocities.v_th_e;
    
    derived_profiles.ln_lambda_ee = collisions.ln_lambda_ee;
    derived_profiles.ln_lambda_ei = collisions.ln_lambda_ei;
    derived_profiles.ln_lambda_ii = collisions.ln_lambda_ii;
    
    % These fields will now work because collisions struct is populated correctly
    derived_profiles.tau_e_calc = collisions.tau_e_calc;
    derived_profiles.tau_e_15   = collisions.tau_e_15;
    derived_profiles.tau_e_17   = collisions.tau_e_17;
    derived_profiles.tau_i_calc = collisions.tau_i_calc;
    derived_profiles.tau_i_15   = collisions.tau_i_15;
    derived_profiles.tau_i_17   = collisions.tau_i_17;
    
    derived_profiles.nu_i_calc  = collisions.nu_i_calc;
    derived_profiles.nu_ei_calc = collisions.nu_ei_calc;
    derived_profiles.nu_i_15    = collisions.nu_i_15;
    derived_profiles.nu_ei_15   = collisions.nu_ei_15;
    
    derived_profiles.nu_star_i = collisionality.nu_star_i;
    derived_profiles.nu_star_e = collisionality.nu_star_e;
    derived_profiles.R_over_Ln = gradients.R_over_Ln;

    fprintf('Derived profiles calculated successfully.\n\n');
end