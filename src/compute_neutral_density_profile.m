function [neutral_profile] = compute_neutral_density_profile(r_fine, constants)
%COMPUTE_NEUTRAL_DENSITY_PROFILE Calculates the neutral hydrogen density profile.
%   This function implements the amplitude version of a Gaussian peak function
%   to model the neutral particle density (n_H0) profile, as described in
%   the thesis (Eq. 3.19). The parameters for the model are sourced from the
%   'constants' structure.
%
%   Syntax:
%       neutral_profile = compute_neutral_density_profile(r_fine, constants)
%
%   Inputs:
%       r_fine    - The high-resolution radial grid [m].
%       constants - Structure containing model parameters under constants.models.neutrals.
%
%   Output:
%       neutral_profile - A structure containing the n_H0 profile.

fprintf('Calculating neutral particle density profile...\n');

% --- Extract model parameters for clarity ---
n_H0_c    = constants.models.neutrals.n_H0_centre;
n_H0_max  = constants.models.neutrals.n_H0_max_amp;
r_max     = constants.models.neutrals.r_max_H_alpha;
width     = constants.models.neutrals.width;

% --- Calculate the profile ---
% This implements the Gaussian peak function from the original monolithic script.
n_H0_profile = n_H0_c + n_H0_max * exp(-0.5 * ((r_fine - r_max) / width).^2);

% --- Package Results into Output Structure ---
neutral_profile.r_fine = r_fine;
neutral_profile.n_H0 = n_H0_profile;

fprintf('Neutral density profile calculated successfully.\n\n');

end