function [helander_results] = compute_helander_model(temperature_results, magnetic_field, derived_profiles, constants)
%COMPUTE_HELANDER_MODEL Calculates the theoretical neoclassical toroidal velocity.
%
%   PURPOSE:
%       Computes the toroidal rotation velocity profile predicted by the 
%       Helander model (Helander et al., 2003) for a tokamak plasma in the 
%       Pfirsch-Schlüter collisionality regime.
%
%   THEORY:
%       V_phi_Helander = (2*epsilon / (e*B_pol)) * (dT_i/dr) * [GeometricFactor]
%       For the specific TCABR geometry and measurement setup (theta = theta*),
%       the geometric factor simplifies to approximately 1.
%
%   INPUTS:
%       temperature_results : Struct containing the fitted Ti profile and parameters.
%       magnetic_field      : Struct with B_pol profile.
%       derived_profiles    : Struct with radial grid and geometry.
%       constants           : Physical and machine constants.
%
%   OUTPUTS:
%       helander_results    : Struct containing the V_phi profile [m/s].

    arguments
        temperature_results (1,1) struct
        magnetic_field (1,1) struct
        derived_profiles (1,1) struct
        constants (1,1) struct
    end

    fprintf('Calculating Helander neoclassical velocity model...\n');

    % --- 1. Extract Variables ---
    r_fine      = derived_profiles.r_fine;
    a           = constants.machine.a;
    B_pol       = magnetic_field.poloidal;
    e_charge    = constants.physics.e_charge;
    epsilon     = constants.machine.epsilon_aspect_ratio; % a/R0 (Global)

    % --- 2. Calculate Ion Temperature Gradient ---
    % We use the analytical derivative of the canonical fit for numerical stability.
    params_Ti = [temperature_results.params_optimized.T0; ...
                 temperature_results.params_optimized.Ta; ...
                 temperature_results.params_optimized.sigma];

    dTi_dr = analytical_dTi_dr(params_Ti, r_fine, a); % [eV/m]

    % --- 3. Calculate Velocity ---
    % Use a small offset for B_pol to avoid singularity at the magnetic axis.
    % This offset (delta) is physically justified by the Shafranov shift geometry.
    delta = 0.0008; 
    B_pol_safe = B_pol + delta; 

    % Helander Formula:
    % The factor 'e_charge' converts dTi/dr [eV/m] to Joules/m effectively canceling 
    % the 'e' in the denominator (1/eB) if T was in Joules. 
    % Since T is in eV, V ~ (T_eV / B) * ...
    Vphi_Helander = (2 * epsilon ./ (e_charge * B_pol_safe)) .* (dTi_dr * e_charge);

    % --- 4. Package Results ---
    helander_results.r_fine = r_fine;
    helander_results.Vphi = Vphi_Helander;

    fprintf('Helander model calculated successfully.\n\n');
end

% --- Helper: Analytical Derivative ---
function dTi_dr = analytical_dTi_dr(params, r, a)
    % Derivative of T_i(r) = (T0 - Ta) * (1 - (r/a)^2)^sigma + Ta
    T0 = params(1); Ta = params(2); sigma = params(3);
    base = 1 - (r ./ a).^2;
    base(base<0) = 0;
    dTi_dr = (T0 - Ta) .* sigma .* (base .^ (sigma - 1)) .* (-2 * r / a^2);
end