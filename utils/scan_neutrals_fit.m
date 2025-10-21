function results = scan_neutrals_fit(velocity_results, temperature_results, collision_profiles, ...
    scaling_variants, scaling_reference, constants, r_fine, amps, widths, opts)
% SCAN_NEUTRALS_FIT  Grid scan of neutral Gaussian parameters (amplitude x width).
%   For each (amp,width) the routine:
%     - constructs n_H0(r) = n0_centre + amp * exp(-(r-r0).^2/(2*w^2))
%     - computes nu_iH0 via compute_nu_iH0_from_nH0
%     - recomputes experimental chi_eff via compute_diffusivity_profile
%     - computes RMSE against the chosen theoretical variant (scaling_reference)
%   The function returns an RMSE map and the best-fit parameters and results.
%
% Minimal required inputs:
%   velocity_results, temperature_results, collision_profiles, scaling_variants,
%   scaling_reference (string), constants, r_fine, amps (vector), widths (vector)
%
% Optional opts fields:
%   opts.centre (radial position of H-alpha peak) default constants.models.neutrals.r_max_H_alpha
%   opts.n0_centre default constants.models.neutrals.n_H0_centre
%   opts.plot_results logical (default true)
%   opts.verbose logical (default true)

if nargin < 11 || isempty(amps), amps = linspace(5e15,2.5e16,12); end
if nargin < 12 || isempty(widths), widths = linspace(0.01,0.04,10); end
if nargin < 13, opts = struct(); end
if ~isfield(opts,'centre'), opts.centre = constants.models.neutrals.r_max_H_alpha; end
if ~isfield(opts,'n0_centre'), opts.n0_centre = constants.models.neutrals.n_H0_centre; end
if ~isfield(opts,'plot_results'), opts.plot_results = true; end
if ~isfield(opts,'verbose'), opts.verbose = true; end

nA = numel(amps); nW = numel(widths);
RMSE = nan(nW, nA);

% Validate scaling_reference
if ~isfield(scaling_variants, scaling_reference)
    error('scan_neutrals_fit: scaling_variants does not contain %s', scaling_reference);
end
theory_vec = scaling_variants.(scaling_reference)(:);

best.rmse = inf;
best.amp = NaN;
best.width = NaN;
best.diff = [];

for iw = 1:nW
    w = widths(iw);
    if opts.verbose, fprintf('scan row %d/%d (width=%.4f m)\n', iw, nW, w); end
    for ia = 1:nA
        amp = amps(ia);
        % construct Gaussian neutral profile
        nH0 = opts.n0_centre + amp .* exp( - (r_fine - opts.centre).^2 ./ (2 * w.^2) );
        % compute nu_iH0 via wrapper (handles bootstrap/means)
        nu_struct = compute_nu_iH0_from_nH0(nH0, temperature_results, constants, r_fine);
        % create local collision_profiles copy and replace nu_iH0
        coll_local = collision_profiles;
        coll_local.nu_iH0 = nu_struct;
        % recompute experimental diffusivity (uses compute_diffusivity_profile)
        diff_local = compute_diffusivity_profile(velocity_results, coll_local, r_fine);
        exp_vec = diff_local.profile_avg(:);
        % choose valid indices
        valid = isfinite(exp_vec) & isfinite(theory_vec) & (theory_vec >= 0);
        if sum(valid) < round(0.5 * numel(valid))
            RMSE(iw,ia) = NaN;
            continue;
        end
        rmse_val = sqrt(mean((exp_vec(valid) - theory_vec(valid)).^2));
        RMSE(iw,ia) = rmse_val;
        if rmse_val < best.rmse
            best.rmse = rmse_val;
            best.amp = amp;
            best.width = w;
            best.diff = diff_local;
            best.nH0 = nH0;
            best.nu = nu_struct;
        end
    end
end

% Package results
results.RMSE = RMSE;
results.amps = amps;
results.widths = widths;
results.best = best;
results.r_fine = r_fine;
results.scaling_reference = scaling_reference;
results.scaling_variants = scaling_variants;

% Plotting
if opts.plot_results
    figure('Name','Neutral scan RMSE','Color','w','Position',[200 200 900 500]);
    imagesc(amps, widths, RMSE);
    set(gca,'YDir','normal');
    colorbar;
    xlabel('n\_H0\_amp [m^{-3}]');
    ylabel('width [m]');
    title(sprintf('RMSE vs neutral params (ref: %s)', scaling_reference));
    hold on;
    plot(best.amp, best.width, 'wx', 'MarkerSize', 12, 'LineWidth', 2);
    hold off;

    % Overlay best-fit χ_eff vs theory
    figure('Name','Best fit overlay','Color','w','Position',[220 220 900 520]);
    % Use plot_chi_eff_comparison if available; otherwise plot directly
    try
        plot_chi_eff_comparison(best.diff, scaling_variants, scaling_variants, constants, r_fine);
    catch
        % Fallback: simple overlay
        r_norm = r_fine / constants.machine.a;
        plot(r_norm, best.diff.profile_avg, 'k-', 'LineWidth', 2, 'DisplayName', 'χ_{φ,eff} (best fit, exp)');
        hold on;
        plot(r_norm, scaling_variants.chi_eff_Rlocal_profile, 'b--', 'LineWidth', 1.6, 'DisplayName', 'χ_{φ,eff} (theory, R_{local})');
        legend('Location','northwest');
        xlabel('r/a'); ylabel('\chi [m^2/s]'); title('Best fit experimental vs theoretical');
        hold off;
    end
end

end
