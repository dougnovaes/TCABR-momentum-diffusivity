function variants = compute_chi_eff_variants(scaling, velocity_results, derived, constants, r_fine)
% Compute several variants of chi_eff from the scaling-law outputs.
% variants: struct with fields .r_norm and fields for each variant (vector)

R0 = constants.machine.R0;
r_norm = r_fine / constants.machine.a;
r_local = R0 * (1 + r_fine / R0);

% base profiles
chi_phi = scaling.chi_phi_solo;        % n x 1
vpinch_solo = scaling.v_pinch_solo;   % n x 1

% compute R/LVphi profile (use velocity poly fit)
Vphi = velocity_results.poly_fit_avg(:);
dVdr = gradient(Vphi, r_fine);
R_over_LVphi_profile = - R0 ./ Vphi .* dVdr;   % may contain NaN where Vphi ~ 0

% evaluation radii indices to test
idx_mid = findClosestIndex(r_fine, constants.machine.a * 0.5);
idx_core = findClosestIndex(r_fine, constants.machine.a * 0.2);
idx_edge = findClosestIndex(r_fine, constants.machine.a * 0.8);

R_over_LVphi_mid = R_over_LVphi_profile(idx_mid);
R_over_Ln_mid = derived.R_over_Ln(idx_mid);

% helper to build chi_eff given R_choice and R_over_LVphi_const
make_chi_eff = @(chi, vpinch, R_choice, R_over_LVphi_const) ...
    chi .* (1 + (R_choice .* vpinch ./ chi) ./ R_over_LVphi_const);

% variants
variants.r_fine = r_fine;
variants.r_norm = r_norm;

% 1) use global R0, R/L at mid
variants.chi_eff_R0_mid = make_chi_eff(chi_phi, vpinch_solo, R0, R_over_LVphi_mid);

% 2) use local R(r) but R/L at mid (hybrid)
variants.chi_eff_Rlocal_mid = make_chi_eff(chi_phi, vpinch_solo, r_local, R_over_LVphi_mid);

% 3) use R0 and R/L evaluated locally (profile form -> vector)
variants.chi_eff_R0_profile = chi_phi .* (1 + (R0 .* vpinch_solo ./ chi_phi) ./ (R_over_LVphi_profile + eps) );

% 4) use local R and R/L profile
variants.chi_eff_Rlocal_profile = chi_phi .* (1 + (r_local .* vpinch_solo ./ chi_phi) ./ (R_over_LVphi_profile + eps) );

% 5) Peeters-style mid-point as in thesis (R/Ln at mid)
variants.chi_eff_peeters_mid = make_chi_eff(chi_phi, vpinch_solo, R0, R_over_Ln_mid); % alternative use

% store the constants used
variants.R_over_LVphi_mid = R_over_LVphi_mid;
variants.R_over_Ln_mid = R_over_Ln_mid;
variants.idx_mid = idx_mid;

end

%% small helper
function idx = findClosestIndex(x, val)
[~, idx] = min(abs(x - val));
end
