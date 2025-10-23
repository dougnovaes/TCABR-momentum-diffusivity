function [neutral_profile] = compute_neutral_density_profile(r_fine, constants)
%COMPUTE_NEUTRAL_DENSITY_PROFILE Calculates the neutral hydrogen density profile.
%   This function implements a Gaussian peak model for the neutral particle
%   density (n_H0) profile, as described in the thesis (Eq. 3.19). The
%   parameters for the model are sourced from the 'constants' structure.
%
%   Model: n_H0(r) = n_H0_c + n_H0_max * exp(-0.5 * ((r - r_max) / width)^2)
%
%   Syntax:
%       neutral_profile = compute_neutral_density_profile(r_fine, constants)

    arguments
        r_fine (:,1) {mustBeNumeric, mustBeReal, mustBeFinite}
        constants (1,1) struct {mustContainFields(constants, {'models'})}
    end

    fprintf('Calculating neutral particle density profile...\n');

    % --- 1. Extract Model Parameters ---
    params = constants.models.neutrals;
    mustContainFields(params, {'n_H0_centre', 'n_H0_max_amp', 'r_max_H_alpha', 'width'});
    
    n_H0_c = params.n_H0_centre;
    n_H0_max = params.n_H0_max_amp;
    r_max = params.r_max_H_alpha;
    width = params.width;

    % --- 2. Calculate the Profile ---
    n_H0_profile = n_H0_c + n_H0_max * exp(-0.5 * ((r_fine - r_max) / width).^2);

    % --- 3. Package Results ---
    neutral_profile.r_fine = r_fine;
    neutral_profile.n_H0 = n_H0_profile;

    fprintf('Neutral density profile calculated successfully.\n\n');
end