function [derived_profiles] = compute_derived_profiles(exp_data, ...
    temperature_results, magnetic_field, constants, r_fine)
%COMPUTE_DERIVED_PROFILES Calculates local and global plasma physics parameters.
%
%   PURPOSE:
%       This function serves as the central physics engine for the analysis pipeline.
%       It computes spatially resolved profiles (local) and volume-averaged 
%       quantities (global) necessary for transport benchmarking.
%
%   METHODOLOGY:
%       1. Local Profiles: Calculates safety factor (q), magnetic shear (s), 
%          thermal velocities, and Coulomb logarithms using standard Spitzer-Braginskii 
%          formulations (Wesson).
%       2. Collisionality: 
%          - Local: Calculates radial profiles of \nu* (Wesson/Bounce) for diagnostics.
%          - Global: Calculates a representative discharge collisionality \nu*_{global} 
%            using Volume-Averaged Density <ne> and Temperature <Te>, following 
%            the definition by Maslov et al. (2009). This scalar is critical for 
%            scaling laws (Solomon) to avoid edge singularities.
%
%   INPUTS:
%       exp_data            : Struct containing raw experimental profiles (ne, ni).
%       temperature_results : Struct with fitted Ti profile and parameters.
%       magnetic_field      : Struct with calculated B-field (toroidal/poloidal).
%       constants           : Project constants (Zeff, machine geometry, etc.).
%       r_fine              : High-resolution radial grid vector [m].
%
%   OUTPUTS:
%       derived_profiles    : Struct containing:
%                             .collisionality (local and global .global_params)
%                             .gradients (R/Ln)
%                             .q_profile, .s_profile
%                             .v_th_i, .v_th_e
%                             .collisions (frequencies and times)
%
%   REFERENCES:
%       - Wesson, "Tokamaks", 3rd Ed. (Collisions, q-profile).
%       - Maslov et al., Nucl. Fusion 49 (2009) (Global effective collisionality).
%       - Solomon et al., Phys. Plasmas 17 (2010) (Scaling law context).
% 
%   UPDATED: Now includes the strict Maslov criterion for density peaking (linear fit).

    arguments
        exp_data (1,1) struct
        temperature_results (1,1) struct
        magnetic_field (1,1) struct
        constants (1,1) struct
        r_fine (:,1) {mustBeNumeric, mustBeReal, mustBeFinite}
    end

    fprintf('Calculating derived physics profiles...\n');

    % --- 1. Extract Core Variables ---
    B_tor   = magnetic_field.toroidal;
    B_pol   = magnetic_field.poloidal;
    a       = constants.machine.a;
    R0      = constants.machine.R0;
    me      = constants.physics.me;
    mp      = constants.physics.mp;
    e       = constants.physics.e_charge;
    Zeff    = constants.plasma.Zeff; 
    
    ne_fine = exp_data.ne_raw;
    ni_fine = exp_data.ni_raw;
    Ti_mean = temperature_results.profile_avg; 
    Te_mean = (509 - Ti_mean(end)) .* (1 - (r_fine / a).^2).^3.3 + Ti_mean(end); 

    % --- 2. Geometry (q, s) ---
    q_and_shear = struct();
    B_pol_safe = B_pol + 1e-9; 
    fq = exp(0.21303 .* r_fine ./ a); 
    q_profile = (r_fine ./ R0) .* (B_tor ./ B_pol_safe) ./ sqrt(1 - (r_fine ./ R0).^2) .* fq;
    q_profile(1) = q_profile(2); 
    s_profile_local = (r_fine ./ q_profile) .* gradient(q_profile, r_fine);
    s_profile_local(1) = 0; 
    q_and_shear.q_profile = q_profile;
    q_and_shear.s_profile = s_profile_local;

    % --- 3. Thermal Velocities ---
    thermal_velocities = struct();
    thermal_velocities.v_th_i = sqrt(2 * e * Ti_mean / mp); 
    thermal_velocities.v_th_e = sqrt(2 * e * Te_mean / me); 

    % --- 4. Collisions (Local) ---
    collisions = struct();
    collisions.ln_lambda_ee = 14.9 - 0.5 * log(ne_fine./1e20) + log(Te_mean./1000);
    collisions.ln_lambda_ei = 15.2 - 0.5 * log(ne_fine ./ 1e20) + log(Te_mean ./ 1000);
    collisions.ln_lambda_ii = 17.3 - 0.5 * log(ni_fine ./ 1e20) + 1.5 * log(Ti_mean ./ 1000);
    collisions.tau_e_calc = 1.09e16 .* (Te_mean./1000).^(3/2) ./ (ni_fine .* collisions.ln_lambda_ei);
    collisions.tau_i_calc = 6.6e17 * sqrt(constants.plasma.m_ion_amu) .* (Ti_mean./1000).^(3/2) ./ (ni_fine .* collisions.ln_lambda_ii);
    % Fixed lambda aliases
    collisions.tau_e_15 = 1.09e16 .* (Te_mean./1000).^(3/2) ./ (ni_fine .* 15);
    collisions.tau_i_15 = 6.6e17 * sqrt(constants.plasma.m_ion_amu) .* (Ti_mean./1000).^(3/2) ./ (ni_fine .* 15);
    collisions.tau_e_17 = 1.09e16 .* (Te_mean./1000).^(3/2) ./ (ni_fine .* 17);
    collisions.tau_i_17 = 6.6e17 * sqrt(constants.plasma.m_ion_amu) .* (Ti_mean./1000).^(3/2) ./ (ni_fine .* 17);
    collisions.nu_ei_15 = 1 ./ collisions.tau_e_15;
    collisions.nu_i_15  = 1 ./ collisions.tau_i_15;
    collisions.nu_ei_calc = 1 ./ collisions.tau_e_calc;
    collisions.nu_i_calc  = 1 ./ collisions.tau_i_calc;

    % --- 5. Collisionality (Local) ---
    collisionality = struct();
    R_coord = R0 * (1 + r_fine / R0);
    epsilon_global = constants.machine.epsilon_aspect_ratio;
    omega_bounce_e = sqrt(epsilon_global) .* thermal_velocities.v_th_e ./ (q_profile .* R_coord);
    omega_bounce_i = sqrt(epsilon_global) .* thermal_velocities.v_th_i ./ (q_profile .* R_coord);
    collisionality.nu_star_e_wesson = collisions.nu_ei_15 ./ (epsilon_global .* omega_bounce_e);
    collisionality.nu_star_i_wesson = collisions.nu_i_15 ./ (epsilon_global .* omega_bounce_i);
    c_s = sqrt(e * Te_mean ./ mp); 
    omega_De = 2 * sqrt(0.1) .* c_s ./ R_coord;
    collisionality.nu_star_e_solomon_local = collisions.nu_ei_15 ./ omega_De;

    % =========================================================================
    % 6. GLOBAL PARAMETERS (MASLOV/SOLOMON DEFINITION)
    % =========================================================================
    fprintf('... calculating global parameters (Restricted Volume Average)\n');
    
    r_norm = r_fine / a;
    
    % A. Define Integration Region (Confinement Zone)
    % Rationale: To characterize the global thermodynamic regime relevant for 
    % transport scaling, we calculate averages over the confinement zone
    % (0.2 <= r/a <= 0.8). This avoids:
    %   1. Core sawteeth/MHD activity (r/a < 0.2)
    %   2. Edge noise and low-density/cold plasma (r/a > 0.8) which would
    %      artificially lower the averages and inflate collisionality.
    mask_zone = (r_norm >= 0.2) & (r_norm <= 0.8); 
    
    % Extract arrays for integration
    r_zone  = r_fine(mask_zone);
    Te_zone = Te_mean(mask_zone);
    ne_zone = ne_fine(mask_zone);
    
    % B. Calculate Zonal Volume Averages
    % Formula: <Y> = Integral(Y * r dr) / Integral(r dr)
    % The weighting factor 'r' accounts for the cylindrical volume element dV.
    integrate_zone = @(y) trapz(r_zone, y .* r_zone);
    vol_norm_zone  = integrate_zone(ones(size(r_zone)));
    
    Te_vol_avg = integrate_zone(Te_zone) / vol_norm_zone;
    ne_vol_avg = integrate_zone(ne_zone) / vol_norm_zone;
    
    % C. Global Collisionality (Maslov Definition)
    % Uses the zonal averages calculated above.
    % Formula: nu_eff = 1e-14 * Zeff * R * <ne> / <Te>^2
    nu_star_global_maslov = 1e-14 * Zeff * R0 * ne_vol_avg / (Te_vol_avg^2);
    
    % D. Global Density Gradient (Maslov Linear Fit)
    % Maslov (2009) defines the characteristic R/Ln as the slope of ln(ne)
    % fitted over the range 0.2 - 0.8.
    mask_grad = (r_norm >= 0.2) & (r_norm <= 0.8);
    
    if sum(mask_grad) > 5
        % Linear fit: ln(n) = - (1/Ln)*r + C  => Slope P(1) = -1/Ln
        P = polyfit(r_fine(mask_grad), log(ne_fine(mask_grad)), 1);
        Ln_inv_maslov = -P(1);
        R_over_Ln_maslov = R0 * Ln_inv_maslov;
    else
        warning('TCABR:Physics', 'Not enough points for Maslov R/Ln fit. Using mid-radius local gradient.');
        [~, idx_mid] = min(abs(r_fine - a*0.5));
        R_over_Ln_maslov = R0 * (-gradient(log(ne_fine), r_fine));
        R_over_Ln_maslov = R_over_Ln_maslov(idx_mid);
    end
    
    % --- Logging & Storage ---
    fprintf('      Zone (0.2-0.8): <Te>=%.2f eV, <ne>=%.2e m^-3\n', Te_vol_avg, ne_vol_avg);
    fprintf('      Global nu* (Maslov): %.4f\n', nu_star_global_maslov);
    fprintf('      Global R/Ln (Maslov): %.4f\n', R_over_Ln_maslov);

    % Store Global Params in output structure
    collisionality.global.Te_avg = Te_vol_avg;
    collisionality.global.ne_avg = ne_vol_avg;
    collisionality.global.nu_star_e_solomon = nu_star_global_maslov;

    % Store Global Params
    collisionality.global.Te_avg = Te_vol_avg;
    collisionality.global.ne_avg = ne_vol_avg;
    collisionality.global.nu_star_e_solomon = nu_star_global_maslov;

    % % --- 6. GLOBAL PARAMETERS (Maslov Definition) ---
    % fprintf('... calculating global parameters (Maslov Criterion)\n');
    % 
    % % A. Volume Averages (Integral de 0 a a)
    % % The weighting "r" in the integral (y .* r_fine) causes the outer 
    % % volume to have much more weight than the core.
    % integrate_vol = @(y) trapz(r_fine, y .* r_fine);
    % vol_norm = integrate_vol(ones(size(r_fine)));
    % Te_vol_avg = integrate_vol(Te_mean) / vol_norm;
    % ne_vol_avg = integrate_vol(ne_fine) / vol_norm;
    % 
    % % B. Global Collisionality
    % nu_star_global_maslov = 1e-14 * Zeff * R0 * ne_vol_avg / (Te_vol_avg^2);
    % 
    % % C. Global Density Gradient (Linear Fit 0.2-0.8)
    % r_norm = r_fine / a;
    % mask_fit = (r_norm >= 0.2) & (r_norm <= 0.8);
    % if sum(mask_fit) > 5
    %     % ln(n) = - (1/Ln)*r + C  => Slope P(1) = -1/Ln
    %     P = polyfit(r_fine(mask_fit), log(ne_fine(mask_fit)), 1);
    %     Ln_inv_maslov = -P(1);
    %     R_over_Ln_maslov = R0 * Ln_inv_maslov;
    % else
    %     warning('Not enough points for Maslov R/Ln fit. Using local mid value.');
    %     [~, idx_mid] = min(abs(r_fine - a*0.5));
    %     R_over_Ln_maslov = R0 * (-gradient(log(ne_fine), r_fine));
    %     R_over_Ln_maslov = R_over_Ln_maslov(idx_mid);
    % end
    % 
    % fprintf('      <Te>: %.2f eV, <ne>: %.2e m^-3\n', Te_vol_avg, ne_vol_avg);
    % fprintf('      Global nu*: %.4f\n', nu_star_global_maslov);
    % fprintf('      Global R/Ln (Maslov): %.4f\n', R_over_Ln_maslov);

    % Store Global Params
    collisionality.global.Te_avg = Te_vol_avg;
    collisionality.global.ne_avg = ne_vol_avg;
    collisionality.global.nu_star_e_solomon = nu_star_global_maslov;
    
    % Gradients Structure
    gradients.R_over_Ln = R0 .* -gradient(log(ne_fine), r_fine); % Local (for ref)
    gradients.R_over_Ln_maslov = R_over_Ln_maslov;             % Global Scalar

    % Aliases
    collisionality.nu_star_e = collisionality.nu_star_e_wesson;
    collisionality.nu_star_i = collisionality.nu_star_i_wesson;
    collisionality.nu_star_e_solomon = collisionality.nu_star_e_solomon_local;

    % --- 7. Package Results ---
    derived_profiles.r_fine     = r_fine;
    derived_profiles.q_profile  = q_profile;
    derived_profiles.s_profile  = s_profile_local;
    derived_profiles.Te_profile = Te_mean;
    derived_profiles.v_th_i     = thermal_velocities.v_th_i;
    derived_profiles.v_th_e     = thermal_velocities.v_th_e;
    derived_profiles.ln_lambda_ee = collisions.ln_lambda_ee;
    derived_profiles.ln_lambda_ei = collisions.ln_lambda_ei;
    derived_profiles.ln_lambda_ii = collisions.ln_lambda_ii;
    derived_profiles.tau_e_calc = collisions.tau_e_calc;
    derived_profiles.tau_e_15 = collisions.tau_e_15;
    derived_profiles.tau_e_17 = collisions.tau_e_17;
    derived_profiles.tau_i_calc = collisions.tau_i_calc;
    derived_profiles.tau_i_15 = collisions.tau_i_15;
    derived_profiles.tau_i_17 = collisions.tau_i_17;
    derived_profiles.nu_i_calc  = collisions.nu_i_calc;
    derived_profiles.nu_ei_calc = collisions.nu_ei_calc;
    derived_profiles.nu_i_15    = collisions.nu_i_15;
    derived_profiles.nu_ei_15   = collisions.nu_ei_15;
    derived_profiles.thermal_velocities = thermal_velocities;
    derived_profiles.collisions = collisions;
    derived_profiles.collisionality = collisionality;
    derived_profiles.nu_star_i = collisionality.nu_star_i;
    derived_profiles.nu_star_e = collisionality.nu_star_e;
    derived_profiles.nu_star_e_solomon = collisionality.nu_star_e_solomon;
    derived_profiles.global_params = collisionality.global;
    derived_profiles.gradients = gradients;
    derived_profiles.R_over_Ln = gradients.R_over_Ln; % Keep local profile for compatibility if needed

    fprintf('Derived profiles calculated successfully.\n\n');
end



% 
%     arguments
%         exp_data (1,1) struct
%         temperature_results (1,1) struct
%         magnetic_field (1,1) struct
%         constants (1,1) struct
%         r_fine (:,1) {mustBeNumeric, mustBeReal, mustBeFinite}
%     end
% 
%     fprintf('Calculating derived physics profiles...\n');
% 
%     % =========================================================================
%     % 1. EXTRACT CORE VARIABLES
%     % =========================================================================
%     B_tor   = magnetic_field.toroidal;
%     B_pol   = magnetic_field.poloidal;
%     a       = constants.machine.a;
%     R0      = constants.machine.R0;
%     me      = constants.physics.me;
%     mp      = constants.physics.mp;
%     e       = constants.physics.e_charge;
%     Zeff    = constants.plasma.Zeff; 
% 
%     ne_fine = exp_data.ne_raw;
%     ni_fine = exp_data.ni_raw;
%     Ti_mean = temperature_results.profile_avg; 
% 
%     % Electron temperature (Canonical shape assumption if experimental Te is noisy)
%     Te_mean = (509 - Ti_mean(end)) .* (1 - (r_fine / a).^2).^3.3 + Ti_mean(end); 
% 
%     % =========================================================================
%     % 2. GEOMETRY: SAFETY FACTOR (q) AND SHEAR (s)
%     % =========================================================================
%     fprintf('... calculating safety factor and magnetic shear\n');
%     q_and_shear = struct();
% 
%     delta = 0.00008; % added to BPol values to avoid the division by ~zero.
%     B_pol_safe = B_pol + delta; 
%     fq = exp(0.21303 .* r_fine ./ a); % fix factor by Severo & Novaes (08/07/2024)
% 
%     q_profile = (r_fine ./ R0) .* (B_tor ./ B_pol_safe) ./ sqrt(1 - (r_fine ./ R0).^2) .* fq;
%     q_profile(1) = q_profile(2); 
% 
%     s_profile_local = (r_fine ./ q_profile) .* gradient(q_profile, r_fine);
%     s_profile_local(1) = 0; 
% 
%     q_and_shear.q_profile = q_profile;
%     q_and_shear.s_profile = s_profile_local;
% 
%     % =========================================================================
%     % 3. THERMAL VELOCITIES
%     % =========================================================================
%     fprintf('... calculating thermal velocities\n');
%     thermal_velocities = struct();
%     thermal_velocities.v_th_i = sqrt(2 * e * Ti_mean / mp); 
%     thermal_velocities.v_th_e = sqrt(2 * e * Te_mean / me); 
% 
%     % =========================================================================
%     % 4. COLLISIONS (LOCAL)
%     % =========================================================================
%     fprintf('... calculating local collision metrics\n');
%     collisions = struct();
% 
%     collisions.ln_lambda_ee = 14.9 - 0.5 * log(ne_fine./1e20) + log(Te_mean./1000);
%     collisions.ln_lambda_ei = 15.2 - 0.5 * log(ne_fine ./ 1e20) + log(Te_mean ./ 1000);
%     collisions.ln_lambda_ii = 17.3 - 0.5 * log(ni_fine ./ 1e20) + 1.5 * log(Ti_mean ./ 1000);
% 
%     collisions.tau_e_calc = 1.09e16 .* (Te_mean./1000).^(3/2) ./ (ni_fine .* collisions.ln_lambda_ei);
%     collisions.tau_i_calc = 6.6e17 * sqrt(constants.plasma.m_ion_amu) .* (Ti_mean./1000).^(3/2) ./ (ni_fine .* collisions.ln_lambda_ii);
% 
%     % Fixed lambda versions for consistency checks
%     collisions.tau_e_15   = 1.09e16 .* (Te_mean./1000).^(3/2) ./ (ni_fine .* 15);
%     collisions.tau_e_17   = 1.09e16 .* (Te_mean./1000).^(3/2) ./ (ni_fine .* 17);
%     collisions.tau_i_15   = 6.6e17 * sqrt(constants.plasma.m_ion_amu) .* (Ti_mean./1000).^(3/2) ./ (ni_fine .* 15);
%     collisions.tau_i_17   = 6.6e17 * sqrt(constants.plasma.m_ion_amu) .* (Ti_mean./1000).^(3/2) ./ (ni_fine .* 17);
% 
%     collisions.nu_ei_calc = 1 ./ collisions.tau_e_calc;
%     collisions.nu_i_calc  = 1 ./ collisions.tau_i_calc;
%     collisions.nu_ei_15   = 1 ./ collisions.tau_e_15;
%     collisions.nu_i_15    = 1 ./ collisions.tau_i_15;
% 
%     % =========================================================================
%     % 5. COLLISIONALITY (LOCAL AND GLOBAL)
%     % =========================================================================
%     fprintf('... calculating collisionality regimes\n');
%     collisionality = struct();
% 
%     R_coord = R0 * (1 + r_fine / R0);
%     epsilon_global = constants.machine.epsilon_aspect_ratio;
% 
%     % --- 5a. Local (Wesson/Bounce) ---
%     omega_bounce_e = sqrt(epsilon_global) .* thermal_velocities.v_th_e ./ (q_and_shear.q_profile .* R_coord);
%     omega_bounce_i = sqrt(epsilon_global) .* thermal_velocities.v_th_i ./ (q_and_shear.q_profile .* R_coord);
% 
%     collisionality.nu_star_e_wesson = collisions.nu_ei_15 ./ (epsilon_global .* omega_bounce_e);
%     collisionality.nu_star_i_wesson = collisions.nu_i_15 ./ (epsilon_global .* omega_bounce_i);
% 
%     % Local Drift (Reference only)
%     c_s = sqrt(e * Te_mean ./ mp); 
%     omega_De = 2 * sqrt(0.1) .* c_s ./ R_coord;
%     collisionality.nu_star_e_solomon_local = collisions.nu_ei_15 ./ omega_De;
% 
%     % --- 5b. Global Representative Collisionality (Maslov/Solomon) ---
%     % Calculates volume-averaged <Te> and <ne> to determine the discharge regime
%     % without edge singularities.
% 
%     % Integration Helper (Cylindrical Shells)
%     integrate_vol = @(y) trapz(r_fine, y .* r_fine);
%     vol_norm = integrate_vol(ones(size(r_fine)));
% 
%     Te_vol_avg = integrate_vol(Te_mean) / vol_norm;
%     ne_vol_avg = integrate_vol(ne_fine) / vol_norm;
% 
%     fprintf('      <Te> (Vol. Avg): %.2f eV\n', Te_vol_avg);
%     fprintf('      <ne> (Vol. Avg): %.2e m^-3\n', ne_vol_avg);
% 
%     % Maslov Engineering Formula: nu_eff = 1e-14 * Zeff * R * <ne> / <Te>^2
%     % <Te> in eV, <ne> in m^-3, R in m.
%     nu_star_global_maslov = 1e-14 * Zeff * R0 * ne_vol_avg / (Te_vol_avg^2);
% 
%     fprintf('      Global Maslov Collisionality (Zeff=%d): %.4f\n', Zeff, nu_star_global_maslov);
% 
%     % Store Global Params
%     collisionality.global.Te_avg = Te_vol_avg;
%     collisionality.global.ne_avg = ne_vol_avg;
%     collisionality.global.nu_star_e_solomon = nu_star_global_maslov;
% 
%     % Aliases for compatibility
%     collisionality.nu_star_e = collisionality.nu_star_e_wesson;
%     collisionality.nu_star_i = collisionality.nu_star_i_wesson;
%     collisionality.nu_star_e_solomon = collisionality.nu_star_e_solomon_local;
% 
%     % =========================================================================
%     % 6. GRADIENTS
%     % =========================================================================
%     gradients = struct();
%     L_ne_inv = -gradient(log(ne_fine), r_fine); 
%     gradients.R_over_Ln = R0 .* L_ne_inv;
% 
%     % =========================================================================
%     % 7. PACKAGE RESULTS
%     % =========================================================================
%     derived_profiles.r_fine     = r_fine;
%     derived_profiles.q_profile  = q_profile;
%     derived_profiles.s_profile  = q_and_shear.s_profile;
%     derived_profiles.Te_profile = Te_mean;
%     derived_profiles.v_th_i     = thermal_velocities.v_th_i;
%     derived_profiles.v_th_e     = thermal_velocities.v_th_e;
% 
%     derived_profiles.ln_lambda_ee = collisions.ln_lambda_ee;
%     derived_profiles.ln_lambda_ei = collisions.ln_lambda_ei;
%     derived_profiles.ln_lambda_ii = collisions.ln_lambda_ii;
%     derived_profiles.tau_e_calc = collisions.tau_e_calc;
%     derived_profiles.tau_e_15 = collisions.tau_e_15;
%     derived_profiles.tau_e_17 = collisions.tau_e_17;
%     derived_profiles.tau_i_calc = collisions.tau_i_calc;
%     derived_profiles.tau_i_15 = collisions.tau_i_15;
%     derived_profiles.tau_i_17 = collisions.tau_i_17;
%     derived_profiles.nu_i_calc  = collisions.nu_i_calc;
%     derived_profiles.nu_ei_calc = collisions.nu_ei_calc;
%     derived_profiles.nu_i_15    = collisions.nu_i_15;
%     derived_profiles.nu_ei_15   = collisions.nu_ei_15;
% 
%     derived_profiles.thermal_velocities = thermal_velocities;
%     derived_profiles.collisions = collisions;
%     derived_profiles.collisionality = collisionality;
% 
%     derived_profiles.nu_star_i = collisionality.nu_star_i;
%     derived_profiles.nu_star_e = collisionality.nu_star_e;
%     derived_profiles.nu_star_e_solomon = collisionality.nu_star_e_solomon;
% 
%     derived_profiles.global_params = collisionality.global;
%     derived_profiles.R_over_Ln  = gradients.R_over_Ln;
% 
%     fprintf('Derived profiles calculated successfully.\n\n');
% end