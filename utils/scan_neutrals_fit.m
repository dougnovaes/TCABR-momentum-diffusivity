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
    % default: one position at constants.models.neutrals.r_max_H_alpha
    r0s = constants.models.neutrals.r_max_H_alpha;
end
if nargin < 9 || isempty(widths), widths = linspace(0.01, 0.04, 10); end
if nargin < 8 || isempty(amps), amps = linspace(5e15, 3e16, 12); end

% Options defaults
if ~isfield(opts,'n0_centre'), opts.n0_centre = constants.models.neutrals.n_H0_centre; end
if ~isfield(opts,'plot_results'), opts.plot_results = true; end
if ~isfield(opts,'verbose'), opts.verbose = true; end
if ~isfield(opts,'eval_rmin_ra'), opts.eval_rmin_ra = 0.2; end
if ~isfield(opts,'eval_rmax_ra'), opts.eval_rmax_ra = 0.99; end

% Convert to column vectors
r_fine = r_fine(:);
amps = amps(:).';
widths = widths(:);
r0s = r0s(:);

nA = numel(amps);
nW = numel(widths);
nR0 = numel(r0s);

% Validate scaling_reference; accept string or cell-array of strings
if ischar(scaling_reference), scaling_fields = {scaling_reference};
elseif iscell(scaling_reference), scaling_fields = scaling_reference;
else error('scaling_reference must be a string or cell array of strings'); end

for k = 1:numel(scaling_fields)
    if ~isfield(scaling_variants, scaling_fields{k})
        error('scan_neutrals_fit: scaling_variants does not contain field "%s"', scaling_fields{k});
    end
end

% Theory vector (primary) - use first field for RMSE compute
theory_primary = scaling_variants.(scaling_fields{1});
theory_primary = theory_primary(:);

% Prepare evaluation index: restrict to r/a in [eval_rmin_ra, eval_rmax_ra]
a = constants.machine.a;
r_norm = r_fine ./ a;
eval_mask = (r_norm >= opts.eval_rmin_ra) & (r_norm <= opts.eval_rmax_ra);
if sum(eval_mask) < 5
    error('scan_neutrals_fit: evaluation mask too small (need at least 5 points). Adjust eval_rmin_ra/eval_rmax_ra.');
end
eval_idx = find(eval_mask);

% Pre-allocate RMSE cube: r0 x width x amp
RMSE = nan(nR0, nW, nA);

best.rmse = inf;
best.amp = NaN; best.width = NaN; best.r0 = NaN;
best.diff = []; best.nH0 = []; best.nu = [];

% ----------------- Main nested scan loop -----------------
for ir0 = 1:nR0
    r0 = r0s(ir0);
    if opts.verbose
        if isnumeric(r0) && isscalar(r0)
            fprintf('Scanning r0 %d/%d (r0=%.4f m)\n', ir0, nR0, r0);
        else
            disp(['Scanning r0 index ', num2str(ir0), '/', num2str(nR0)]);
        end
    end
    for iw = 1:nW
        w = widths(iw);
        for ia = 1:nA
            amp = amps(ia);
            % construct Gaussian neutral profile
            nH0 = opts.n0_centre + amp .* exp( - (r_fine - r0).^2 ./ (2 * w.^2) );
            % compute nu_iH0 via wrapper
            nu_struct = compute_nu_iH0_from_nH0(nH0, temperature_results, constants, r_fine);
            % create local collision profiles and replace nu_iH0
            coll_local = collision_profiles;
            coll_local.nu_iH0 = nu_struct;
            % recompute experimental diffusivity
            diff_local = compute_diffusivity_profile(velocity_results, coll_local, r_fine);
            exp_vec = diff_local.profile_avg(:);
            % compute RMSE comparing against theory_primary restricted to eval_idx
            tvec = theory_primary(:);
            % select valid indices inside eval_idx excluding NaNs or negatives in theory
            valid = eval_mask & isfinite(exp_vec) & isfinite(tvec) & (tvec >= 0);
            if sum(valid) < round(0.5 * sum(eval_mask))
                RMSE(ir0, iw, ia) = NaN;
                continue;
            end
            % compute RMSE ignoring NaNs
            diff = exp_vec(valid) - tvec(valid);
            rmse_val = sqrt(mean(diff.^2));
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

% ----------------- Package results -----------------
results.RMSE = RMSE;
results.amps = amps;
results.widths = widths;
results.r0s = r0s;
results.best = best;
results.r_fine = r_fine;
results.eval_mask = eval_mask;
results.scaling_reference = scaling_reference;
results.scaling_variants = scaling_variants;

% ----------------- Plotting -----------------
if opts.plot_results
    % 1) RMSE best-slice: show best r0 plane as heatmap (width x amp)
    [~,ir0_best] = min(reshape(nanmin(reshape(RMSE,[],nA),1),[],1)); %#ok - fallback
    % but better use best.r0
    if ~isnan(best.r0)
        ir0_best = find(abs(r0s - best.r0) == min(abs(r0s - best.r0)),1);
    else
        ir0_best = 1;
    end

    figure('Name','Neutral scan RMSE (best r0)','Color','w','Position',[200 200 900 500]);
    rmse_plane = squeeze(RMSE(ir0_best,:,:)); % width x amp
    imagesc(amps, widths, rmse_plane);
    set(gca,'YDir','normal');
    xlabel('n\_H0\_amp [m^{-3}]'); ylabel('width [m]');
    title(sprintf('RMSE vs neutral params (r0 = %.3f m)', r0s(ir0_best)));
    colorbar;
    hold on;
    plot(best.amp, best.width, 'wx', 'MarkerSize', 12, 'LineWidth', 2);
    % annotate chosen eval window and theory primary name and mid values
    % compute and display R/L_n(mid), R/L_Vphi(mid) if available in scaling_variants
    mid_idx = find(abs(r_fine - a*0.5) == min(abs(r_fine - a*0.5)),1);
    text(0.02*max(amps)+min(amps), 0.95*max(widths)+min(widths), ...
        sprintf('R0=%.3f m\nR/L_n(mid)=%.2f', constants.machine.R0, ...
        (isfield(results.scaling_variants,'R_over_Ln_const')*results.scaling_variants.R_over_Ln_const)), ...
        'BackgroundColor','w','EdgeColor','k','FontSize',9);

    hold off;

    % 2) Overlay plot: experimental best vs multiple theoretical variants
    figure('Name','Best fit overlay','Color','w','Position',[220 220 900 520]);
    r_norm = r_fine ./ a;
    % plot experimental best with CI
    plot(r_norm, best.diff.profile_avg, 'k-', 'LineWidth', 2, 'DisplayName', '\chi_{eff} (Exp, best)');
    hold on;
    if isfield(best.diff,'ci_lower') && isfield(best.diff,'ci_upper')
        fill([r_norm; flipud(r_norm)], [best.diff.ci_lower; flipud(best.diff.ci_upper)], ...
            [0.6 0.6 0.6], 'FaceAlpha', 0.25, 'EdgeColor', 'none', 'HandleVisibility','off');
    end

    % Plot each requested theory variant with explicit names
    colors = get(gca,'ColorOrder');
    ncols = size(colors,1);
    for k = 1:numel(scaling_fields)
        fld = scaling_fields{k};
        th = scaling_variants.(fld)(:);
        % Safe plotting: mask where eval_mask true OR plot full line with dashed where invalid
        h = plot(r_norm, th, 'LineWidth', 1.6, 'DisplayName', sprintf('%s (theory)', fld));
        % cycle colours manually to preserve variety
        set(h, 'Color', colors(mod(k-1,ncols)+1,:));
    end

    xlabel('r/a'); ylabel('\chi_{\phi,eff} [m^2/s]');
    title('Best-fit experimental effective diffusivity vs theoretical variants');
    legend('Location','northwest','Interpreter','none');
    xlim([0 1]);
    % show evaluation window as vertical shaded region
    ylims = ylim;
    xpatch = [opts.eval_rmin_ra opts.eval_rmax_ra opts.eval_rmax_ra opts.eval_rmin_ra];
    patch(a*xpatch, [ylims(1) ylims(1) ylims(2) ylims(2)], [0.9 0.9 0.9], 'FaceAlpha', 0.15, 'EdgeColor','none', 'HandleVisibility','off');
    text( (opts.eval_rmin_ra+opts.eval_rmax_ra)/2, ylims(2)*0.95, 'RMSE eval window', 'HorizontalAlignment','center', 'BackgroundColor','w','EdgeColor','k','FontSize',9 );
    hold off;
end

end
