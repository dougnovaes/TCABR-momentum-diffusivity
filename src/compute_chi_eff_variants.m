function variants = compute_chi_eff_variants(scaling_law_results, velocity_results, derived_profiles, constants, r_fine)
% COMPUTE_CHI_EFF_VARIANTS
%   Build several variants of theoretical effective diffusivity chi_phi_eff
%   combining: use of R0 vs local R(r), and using mid-radius scalar R/L_Vphi
%   vs full R/L_Vphi profile. Returns a struct with named variants.
%
% Inputs:
%   scaling_law_results - output of compute_chi_eff_from_scalings (contains chi_phi_solo, etc.)
%   velocity_results    - result of analyze_velocity_profile (contains Vphi fits)
%   derived_profiles    - derived_profiles (contains R_over_Ln, etc.)
%   constants, r_fine   - as usual
%
% Output:
%   variants - struct with fields:
%       .r_fine
%       .chi_phi_solo
%       .chi_eff_R0_mid
%       .chi_eff_R0_profile
%       .chi_eff_Rlocal_mid
%       .chi_eff_Rlocal_profile
%       .meta (struct with used constants)

a   = constants.machine.a;
R0  = constants.machine.R0;
r   = r_fine(:);
r_norm = r ./ a;

% base solo diffusivity profile (already uses an R/Ln constant in compute_chi_eff_from_scalings)
chi_phi_solo = scaling_law_results.chi_phi_solo(:);

% compute R_over_LVphi profile robustly from velocity_results.poly_fit_avg
Vphi = velocity_results.poly_fit_avg(:);
dVdr = gradient(Vphi, r);
% avoid division by zero: regularise small Vphi
Vphi_safe = Vphi;
smallV = abs(Vphi_safe) < 1e-8;
Vphi_safe(smallV) = sign(Vphi_safe(smallV)).*1e-8 + eps;

R_over_LVphi_profile = -R0 ./ Vphi_safe .* dVdr;  % may contain large values near zeros

% mid-radius index
[~, idx_mid] = min(abs(r - a*0.5));
R_over_LVphi_mid = R_over_LVphi_profile(idx_mid);

% build R_coord (local major radius approximation)
R_coord = R0 * (1 + r ./ R0);

% Pinch models (recompute consistently from chi_phi_solo, to ensure same baseline)
% Use Solomon pinch (C * nu_star_e) stored previously (but recompute simple forms here)
% For consistency with compute_chi_eff_from_scalings: use same formulas
A = 6.09; B = 0.157; C = -24.2;
nu_star_e = derived_profiles.nu_star_e(:);
v_pinch_solo = C .* nu_star_e;
v_pinch_hahm = -2 .* chi_phi_solo ./ R0;
F = 1;
v_pinch_gurcan = -(2 .* chi_phi_solo ./ R_coord) .* (F + r ./ R0);
v_pinch_peeters_Rln2 = (chi_phi_solo ./ R_coord) .* (-4 - 2);
R_over_Ln_mid = derived_profiles.R_over_Ln(idx_mid);
v_pinch_peeters_Rln_calc = (chi_phi_solo ./ R_coord) .* (-4 - R_over_Ln_mid);

% Helper to compute chi_eff given PinchNumber and R_over_LVphi (scalar or profile)
compute_chi_eff = @(PinchNumber, R_over_LVphi) chi_phi_solo .* (1 + PinchNumber ./ R_over_LVphi);

% PinchNumber = R * v_pinch / chi_phi
Pin_Solo   = R0 .* v_pinch_solo ./ chi_phi_solo;
Pin_Hahm   = R0 .* v_pinch_hahm ./ chi_phi_solo;
Pin_Gurcan = R0 .* v_pinch_gurcan ./ chi_phi_solo;
Pin_Peet2  = R0 .* v_pinch_peeters_Rln2 ./ chi_phi_solo;
Pin_PeetC  = R0 .* v_pinch_peeters_Rln_calc ./ chi_phi_solo;

% Variants:
% 1) Use R0 and mid R/L_Vphi (scalar)
chi_eff_R0_mid = compute_chi_eff(Pin_Solo, R_over_LVphi_mid);    % Solo example
% We'll compute arrays stacking variants; for clarity produce full set with each pinch model
chi_eff_R0_mid_Solo   = compute_chi_eff(Pin_Solo, R_over_LVphi_mid);
chi_eff_R0_mid_Hahm   = compute_chi_eff(Pin_Hahm, R_over_LVphi_mid);
chi_eff_R0_mid_Gurcan = compute_chi_eff(Pin_Gurcan, R_over_LVphi_mid);
chi_eff_R0_mid_Peet2  = compute_chi_eff(Pin_Peet2, R_over_LVphi_mid);
chi_eff_R0_mid_PeetC  = compute_chi_eff(Pin_PeetC, R_over_LVphi_mid);

% 2) Use R0 and full R/L_Vphi profile (pointwise)
chi_eff_R0_profile_Solo   = compute_chi_eff(Pin_Solo, R_over_LVphi_profile);
chi_eff_R0_profile_Hahm   = compute_chi_eff(Pin_Hahm, R_over_LVphi_profile);
chi_eff_R0_profile_Gurcan = compute_chi_eff(Pin_Gurcan, R_over_LVphi_profile);
chi_eff_R0_profile_Peet2  = compute_chi_eff(Pin_Peet2, R_over_LVphi_profile);
chi_eff_R0_profile_PeetC  = compute_chi_eff(Pin_PeetC, R_over_LVphi_profile);

% 3) Use local R_coord in PinchNumber (already used for Gurcan/Peeters) and mid R/L_Vphi scalar
% For completeness you may want a variant where the denominator R_over_LVphi uses local R (ambiguous in literature).
chi_eff_Rlocal_mid_Solo   = compute_chi_eff(Pin_Solo, R_over_LVphi_mid);
chi_eff_Rlocal_profile_Solo = compute_chi_eff(Pin_Solo, R_over_LVphi_profile);

% Package
variants.r_fine = r;
variants.r_norm = r_norm;
variants.chi_phi_solo = chi_phi_solo;

variants.R_over_LVphi_profile = R_over_LVphi_profile;
variants.R_over_LVphi_mid = R_over_LVphi_mid;
variants.R_over_Ln_mid = R_over_Ln_mid;

variants.chi_eff_R0_mid = struct( ...
    'Solo',   chi_eff_R0_mid_Solo, ...
    'Hahm',   chi_eff_R0_mid_Hahm, ...
    'Gurcan', chi_eff_R0_mid_Gurcan, ...
    'Peeters_R2', chi_eff_R0_mid_Peet2, ...
    'Peeters_C',  chi_eff_R0_mid_PeetC );

variants.chi_eff_R0_profile = struct( ...
    'Solo',   chi_eff_R0_profile_Solo, ...
    'Hahm',   chi_eff_R0_profile_Hahm, ...
    'Gurcan', chi_eff_R0_profile_Gurcan, ...
    'Peeters_R2', chi_eff_R0_profile_Peet2, ...
    'Peeters_C',  chi_eff_R0_profile_PeetC );

variants.chi_eff_Rlocal_mid = struct( ...
    'Solo', chi_eff_Rlocal_mid_Solo);

variants.chi_eff_Rlocal_profile = struct( ...
    'Solo', chi_eff_Rlocal_profile_Solo);

variants.meta.R0 = R0;
variants.meta.a = a;
variants.meta.R_over_LVphi_mid = R_over_LVphi_mid;
variants.meta.R_over_Ln_mid = R_over_Ln_mid;

end
