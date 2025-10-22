function results = scan_neutrals_fit(velocity_results, temperature_results, collision_profiles, ...
    scaling_variants, scaling_reference, constants, r_fine, amps, widths, r0s, opts)
% SCAN_NEUTRALS_FIT  Grid scan of neutral Gaussian parameters (amp × width × r0).
%   Extended version: also scans the Gaussian peak position r0 and evaluates
%   RMSE only on a safe radial window (r/a ∈ [0.2, 0.99]) to avoid divergences
%   where 1/V or gradV/V produce singularities.
%
% USAGE:
%   results = scan_neutrals_fit(velocity_results, temperature_results, collision_profiles, ...
%       scaling_variants, scaling_reference, constants, r_fine, amps, widths, r0s, opts)
%
% REQUIRED INPUTS:
%   velocity_results, temperature_results, collision_profiles
%   scaling_variants  - struct containing theoretical variant vectors (fields)
%   scaling_reference - name of theory field to compare (string) OR cell array of fields to overlay
%   constants, r_fine  - as usual
%   amps, widths, r0s  - vectors of amplitudes [m^-3], widths [m], peak positions [m]
%
% OPTIONAL 'opts' fields:
%   opts.n0_centre    - background neutral density [m^-3] (default from constants)
%   opts.plot_results - logical (default true)
%   opts.verbose      - logical (default true)
%   opts.eval_rmin_ra - minimum r/a for RMSE eval (default 0.2)
%   opts.eval_rmax_ra - maximum r/a for RMSE eval (default 0.99)
%
% OUTPUT:
%   results struct containing RMSE cube (nR0 x nW x nA), best parameters, and figures.

% ----------------- Defaults & input checks -----------------
if nargin < 11, opts = struct(); end
if nargin < 10 || isempty(r0s)
    r0s = constants.models.neutrals.r_max_H_alpha;
end
if nargin < 9 || isempty(widths), widths = linspace(0.01,0.04,10); end
if nargin < 8 || isempty(amps), amps = linspace(5e15,3e16,12); end

% options defaults
if ~isfield(opts,'n0_centre'), opts.n0_centre = constants.models.neutrals.n_H0_centre; end
if ~isfield(opts,'plot_results'), opts.plot_results = true; end
if ~isfield(opts,'verbose'), opts.verbose = true; end
if ~isfield(opts,'eval_rmin_ra'), opts.eval_rmin_ra = 0.2; end
if ~isfield(opts,'eval_rmax_ra'), opts.eval_rmax_ra = 0.99; end

% prepare sizes and preallocate RMSE cube: nR0 x nW x nA
amps = amps(:).'; widths = widths(:); r0s = r0s(:);
nA = numel(amps); nW = numel(widths); nR0 = numel(r0s);
RMSE = nan(nR0, nW, nA);

% choose primary theory vector
if ischar(scaling_reference), scaling_reference = {scaling_reference}; end
if ~isfield(scaling_variants, scaling_reference{1})
    error('scan_neutrals_fit: scaling_variants does not contain %s', scaling_reference{1});
end
theory_primary = scaling_variants.(scaling_reference{1});
theory_primary = theory_primary(:);

% evaluation mask
a = constants.machine.a;
r_norm = r_fine ./ a;
eval_mask = (r_norm >= opts.eval_rmin_ra) & (r_norm <= opts.eval_rmax_ra);
eval_idx = find(eval_mask);
if numel(eval_idx) < 5
    error('scan_neutrals_fit: evaluation mask too small.');
end

% prepare best struct
best.rmse = inf;
best.amp = NaN; best.width = NaN; best.r0 = NaN;
best.diff = []; best.nH0 = []; best.nu = [];

% Main loops (preallocated RMSE prevents growing)
for ir0 = 1:nR0
    r0 = r0s(ir0);
    if opts.verbose
        fprintf('Scanning r0 %d/%d (r0=%.4f m)\n', ir0, nR0, r0);
    end
    for iw = 1:nW
        w = widths(iw);
        for ia = 1:nA
            amp = amps(ia);
            % construct neutral profile
            nH0 = opts.n0_centre + amp .* exp( - (r_fine - r0).^2 ./ (2 * w.^2) );
            % compute nu using wrapper
            nu_struct = compute_nu_iH0_from_nH0(nH0, temperature_results, constants, r_fine);
            coll_local = collision_profiles;
            coll_local.nu_iH0 = nu_struct;
            diff_local = compute_diffusivity_profile(velocity_results, coll_local, r_fine);
            exp_vec = diff_local.profile_avg(:);
            tvec = theory_primary(:);
            valid_mask = eval_mask & isfinite(exp_vec) & isfinite(tvec) & (tvec >= 0);
            if sum(valid_mask) < round(0.5 * sum(eval_mask))
                RMSE(ir0, iw, ia) = NaN;
                continue;
            end
            d = exp_vec(valid_mask) - tvec(valid_mask);
            rmse_val = sqrt(mean(d.^2));
            RMSE(ir0, iw, ia) = rmse_val;
            if rmse_val < best.rmse
                best.rmse = rmse_val;
                best.amp = amp;
                best.width = w;
                best.r0 = r0;
                best.diff = diff_local;
                best.nH0 = nH0;
                best.nu = nu_struct;
            end
        end
    end
end

% package
results.RMSE = RMSE;
results.amps = amps;
results.widths = widths;
results.r0s = r0s;
results.best = best;
results.r_fine = r_fine;
results.eval_mask = eval_mask;
results.scaling_reference = scaling_reference;
results.scaling_variants = scaling_variants;

% plotting (kept similar to previous)
if opts.plot_results
    % pick best r0 plane
    if ~isnan(best.r0)
        [~, ir0_best] = min(abs(r0s - best.r0));
    else
        ir0_best = 1;
    end
    figure('Name','Neutral scan RMSE (best r0)','Color','w','Position',[200 200 900 500]);
    rmse_plane = squeeze(RMSE(ir0_best,:,:));
    imagesc(amps, widths, rmse_plane); set(gca,'YDir','normal'); colorbar;
    xlabel('n\_H0\_amp [m^{-3}]'); ylabel('width [m]');
    title(sprintf('RMSE vs neutral params (r0 = %.3f m)', r0s(ir0_best)));
    hold on;
    plot(best.amp, best.width, 'wx', 'MarkerSize', 12, 'LineWidth', 2);
    hold off;

    % overlay best fit
    figure('Name','Best fit overlay','Color','w','Position',[220 220 900 520]);
    r_norm = r_fine ./ a;
    plot(r_norm, best.diff.profile_avg, 'k-', 'LineWidth', 2, 'DisplayName', '\chi_{eff} (best fit, exp)');
    hold on;
    if isfield(best.diff,'ci_lower') && isfield(best.diff,'ci_upper')
        fill([r_norm; flipud(r_norm)], [best.diff.ci_lower; flipud(best.diff.ci_upper)], [0.85 0.85 0.85], 'FaceAlpha', 0.25, 'EdgeColor','none');
    end
    plot(r_norm, scaling_variants.(scaling_reference{1})(:), 'b--', 'LineWidth', 1.6, 'DisplayName', scaling_reference{1});
    legend('Location','northwest');
    xlabel('r/a'); ylabel('\chi [m^2/s]');
    title('Best fit experimental vs theoretical');
    hold off;
end
end
