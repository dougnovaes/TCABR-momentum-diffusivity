function results = scan_neutrals_fit(velocity_results, temperature_results, ...
    variants, constants, r_fine, amps, widths, r0s, opts)
%SCAN_NEUTRALS_FIT Performs a grid scan of neutral Gaussian parameters.
%   Finds the optimal parameters (amplitude, width, peak position) for the
%   neutral density profile by minimizing the RMSE between the resulting
%   experimental chi_phi_eff and a specified theoretical variant.
%   The scan loops are linearized for efficient parallel execution.

    arguments
        velocity_results (1,1) struct
        temperature_results (1,1) struct
        variants (1,1) struct % Correct: Expects the theoretical_variants struct
        constants (1,1) struct
        r_fine (:,1) double
        amps (1,:) double
        widths (1,:) double
        r0s (1,:) double
        opts (1,1) struct
    end

    % --- 1. Setup ---
    % This line should now work correctly as 'variants' is the correct struct
    safe_mask = variants.meta.safe_mask;
    
    ref_parts = split(opts.reference_theory, '.');
    try
        theory_ref_struct = variants.(ref_parts{1}).(ref_parts{2});
        theory_primary = theory_ref_struct.profile_avg;
    catch ME
        error('Invalid reference theory string: ''%s''. Check format (e.g., ''Solomon.R0_mid'') and ensure variant exists. Details: %s', opts.reference_theory, ME.message);
    end
    
    % --- 2. Linearize Scan Parameters for Efficient Parallelization ---
    [R0_grid, W_grid, A_grid] = ndgrid(r0s, widths, amps);
    params_list = [R0_grid(:), W_grid(:), A_grid(:)];
    num_calcs = size(params_list, 1);
    rmse_vector = nan(num_calcs, 1);
    
    fprintf('Starting grid scan (%d total points) using parallel pool...\n', num_calcs);

    parfor i = 1:num_calcs  
        r0  = params_list(i, 1);
        w   = params_list(i, 2);
        amp = params_list(i, 3);
        
        temp_constants = constants;
        temp_constants.models.neutrals.n_H0_max_amp = amp;
        temp_constants.models.neutrals.width = w;
        temp_constants.models.neutrals.r_max_H_alpha = r0;
        
        neutral_prof = compute_neutral_density_profile(r_fine, temp_constants);
        coll_prof = compute_collision_profiles(temperature_results, neutral_prof, r_fine);
        diff_exp = compute_diffusivity_profile(velocity_results, coll_prof, r_fine);
        
        residuals = diff_exp.profile_avg(safe_mask) - theory_primary(safe_mask);
        rmse_vector(i) = sqrt(mean(residuals.^2));
    end
    
    % --- 3. Find Best Result and Reshape RMSE data ---
    RMSE_cube = reshape(rmse_vector, numel(r0s), numel(widths), numel(amps));
    
    [min_rmse, min_idx] = min(RMSE_cube(:));
     if isnan(min_rmse)
         warning('Scan failed to produce valid RMSE values. Check input ranges or theory reference.');
         best.rmse = inf; best.amp = NaN; best.width = NaN; best.r0 = NaN;
     else
        [idx_r0, idx_w, idx_a] = ind2sub(size(RMSE_cube), min_idx);
        best = struct();
        best.rmse = min_rmse;
        best.amp = amps(idx_a);
        best.width = widths(idx_w);
        best.r0 = r0s(idx_r0);
        
        fprintf('Scan complete. Best RMSE = %.4f\n', best.rmse);
        fprintf('  Best params: Amp=%.2e, Width=%.4f, r0=%.4f\n', best.amp, best.width, best.r0);
     end
     
    % --- 4. Package Results ---
    results.opts = opts;
    results.RMSE_cube = RMSE_cube;
    results.amps = amps;
    results.widths = widths;
    results.r0s = r0s;
    results.best = best;

    % --- 5. Plotting (Optional) ---
    if opts.plot_results && ~isnan(best.rmse)
        fprintf('Generating study-specific plots...\n');
        
        % Plotting logic: Heatmap and Best fit overlay
        try
            a = constants.machine.a; % Get minor radius
            r_norm = r_fine / a;     % Define r_norm for plotting
            % Best r0 plane for heatmap
            [~, ir0_best] = min(abs(r0s - best.r0));
            
            fig_hm = figure('Name','Neutral Scan RMSE (Best r0 Plane)','Color','w','Position',[200 200 900 500], 'WindowStyle', 'docked');
            rmse_plane = squeeze(RMSE_cube(ir0_best,:,:));
            imagesc(amps, widths, rmse_plane); 
            set(gca,'YDir','normal'); colorbar;
            xlabel('Neutral Peak Amplitude, $n_{H0,max}$ [m$^{-3}$]', 'Interpreter', 'latex'); 
            ylabel('Gaussian Width, $w$ [m]', 'Interpreter', 'latex');
            title(sprintf('RMSE vs Neutral Params (at $r_{peak} = %.4f$ m)', r0s(ir0_best)), 'Interpreter', 'latex');
            hold on;
            plot(best.amp, best.width, 'wx', 'MarkerSize', 12, 'LineWidth', 2, 'DisplayName', 'Best Fit');
            hold off;
            set_publication_style(gca);
            
            % Overlay best fit experimental vs theoretical reference
            fig_overlay = figure('Name','Neutral Scan Best Fit Overlay','Color','w','Position',[220 220 900 520], 'WindowStyle', 'docked');
            ax_overlay = gca; hold(ax_overlay, 'on');
            
            % Recompute best experimental profile for plotting CIs
            best_constants = constants;
            best_constants.models.neutrals.n_H0_max_amp = best.amp;
            best_constants.models.neutrals.width = best.width;
            best_constants.models.neutrals.r_max_H_alpha = best.r0;
            best_neutral_prof = compute_neutral_density_profile(r_fine, best_constants);
            best_coll_prof = compute_collision_profiles(temperature_results, best_neutral_prof, r_fine);
            best_diff_exp = compute_diffusivity_profile(velocity_results, best_coll_prof, r_fine);

            plot(ax_overlay, r_norm, best_diff_exp.profile_avg, 'k-', 'LineWidth', 2.5, 'DisplayName', '$\chi_{\phi, \mathrm{eff}}^{\mathrm{(exp)}}$ (Best Fit)');
            fill(ax_overlay, [r_norm; flipud(r_norm)], [best_diff_exp.ci_lower; flipud(best_diff_exp.ci_upper)], ...
                 'k', 'FaceAlpha', 0.15, 'EdgeColor','none', 'DisplayName', '95% CI (Exp.)');
            
            plot(ax_overlay, r_norm, theory_primary, 'b--', 'LineWidth', 2.0, 'DisplayName', ['Theory: ', strrep(opts.reference_theory, '_', '\_')]);
            
            legend(ax_overlay, 'Location','northwest');
            xlabel(ax_overlay, 'Normalised Radius ($r/a$)'); 
            ylabel(ax_overlay, 'Effective Diffusivity, $\chi_{\phi, \mathrm{eff}}$ [m$^2$/s]');
            title(ax_overlay, 'Best Fit Comparison: Experiment vs. Theory');
            xlim(ax_overlay, [0, 1]); ylim(ax_overlay, 'auto');
            grid(ax_overlay, 'on'); box(ax_overlay, 'on');
            set_publication_style(ax_overlay);
            hold(ax_overlay, 'off');

            % Save figures
            scan_plots_dir = fullfile(results_dir, 'studies_plots');
            if ~exist(scan_plots_dir, 'dir'), mkdir(scan_plots_dir); end
            saveas(fig_hm, fullfile(scan_plots_dir, 'neutral_scan_heatmap.png'));
            savefig(fig_hm, fullfile(scan_plots_dir, 'neutral_scan_heatmap.fig'));
            saveas(fig_overlay, fullfile(scan_plots_dir, 'neutral_scan_best_overlay.png'));
            savefig(fig_overlay, fullfile(scan_plots_dir, 'neutral_scan_best_overlay.fig'));
            
        catch ME_plot
            warning(ME_plot.identifier, 'scan_neutrals_fit: failed to generate study plots — %s', ME_plot.message);
        end

    elseif opts.plot_results && isnan(best.rmse)
        warning('Skipping study plots because scan failed to find a valid minimum RMSE.');
    end
end

% function results = scan_neutrals_fit(velocity_results, temperature_results, collision_profiles, ...
%     scaling_variants, scaling_reference, constants, r_fine, amps, widths, r0s, opts)
% % SCAN_NEUTRALS_FIT  Grid scan of neutral Gaussian parameters (amp × width × r0).
% %   Extended version: also scans the Gaussian peak position r0 and evaluates
% %   RMSE only on a safe radial window (r/a ∈ [0.2, 0.99]) to avoid divergences
% %   where 1/V or gradV/V produce singularities.
% %
% % USAGE:
% %   results = scan_neutrals_fit(velocity_results, temperature_results, collision_profiles, ...
% %       scaling_variants, scaling_reference, constants, r_fine, amps, widths, r0s, opts)
% %
% % REQUIRED INPUTS:
% %   velocity_results, temperature_results, collision_profiles
% %   scaling_variants  - struct containing theoretical variant vectors (fields)
% %   scaling_reference - name of theory field to compare (string) OR cell array of fields to overlay
% %   constants, r_fine  - as usual
% %   amps, widths, r0s  - vectors of amplitudes [m^-3], widths [m], peak positions [m]
% %
% % OPTIONAL 'opts' fields:
% %   opts.n0_centre    - background neutral density [m^-3] (default from constants)
% %   opts.plot_results - logical (default true)
% %   opts.verbose      - logical (default true)
% %   opts.eval_rmin_ra - minimum r/a for RMSE eval (default 0.2)
% %   opts.eval_rmax_ra - maximum r/a for RMSE eval (default 0.99)
% %
% % OUTPUT:
% %   results struct containing RMSE cube (nR0 x nW x nA), best parameters, and figures.
% 
% % ----------------- Defaults & input checks -----------------
% if nargin < 11, opts = struct(); end
% if nargin < 10 || isempty(r0s)
%     r0s = constants.models.neutrals.r_max_H_alpha;
% end
% if nargin < 9 || isempty(widths), widths = linspace(0.01,0.04,10); end
% if nargin < 8 || isempty(amps), amps = linspace(5e15,3e16,12); end
% 
% % options defaults
% if ~isfield(opts,'n0_centre'), opts.n0_centre = constants.models.neutrals.n_H0_centre; end
% if ~isfield(opts,'plot_results'), opts.plot_results = true; end
% if ~isfield(opts,'verbose'), opts.verbose = true; end
% if ~isfield(opts,'eval_rmin_ra'), opts.eval_rmin_ra = 0.2; end
% if ~isfield(opts,'eval_rmax_ra'), opts.eval_rmax_ra = 0.99; end
% 
% % prepare sizes and preallocate RMSE cube: nR0 x nW x nA
% amps = amps(:).'; widths = widths(:); r0s = r0s(:);
% nA = numel(amps); nW = numel(widths); nR0 = numel(r0s);
% RMSE = nan(nR0, nW, nA);
% 
% % choose primary theory vector
% if ischar(scaling_reference), scaling_reference = {scaling_reference}; end
% if ~isfield(scaling_variants, scaling_reference{1})
%     error('scan_neutrals_fit: scaling_variants does not contain %s', scaling_reference{1});
% end
% theory_primary = scaling_variants.(scaling_reference{1});
% theory_primary = theory_primary(:);
% 
% % evaluation mask
% a = constants.machine.a;
% r_norm = r_fine ./ a;
% eval_mask = (r_norm >= opts.eval_rmin_ra) & (r_norm <= opts.eval_rmax_ra);
% eval_idx = find(eval_mask);
% if numel(eval_idx) < 5
%     error('scan_neutrals_fit: evaluation mask too small.');
% end
% 
% % prepare best struct
% best.rmse = inf;
% best.amp = NaN; best.width = NaN; best.r0 = NaN;
% best.diff = []; best.nH0 = []; best.nu = [];
% 
% % Main loops (preallocated RMSE prevents growing)
% for ir0 = 1:nR0
%     r0 = r0s(ir0);
%     if opts.verbose
%         fprintf('Scanning r0 %d/%d (r0=%.4f m)\n', ir0, nR0, r0);
%     end
%     for iw = 1:nW
%         w = widths(iw);
%         for ia = 1:nA
%             amp = amps(ia);
%             % construct neutral profile
%             nH0 = opts.n0_centre + amp .* exp( - (r_fine - r0).^2 ./ (2 * w.^2) );
%             % compute nu using wrapper
%             nu_struct = compute_nu_iH0_from_nH0(nH0, temperature_results, constants, r_fine);
%             coll_local = collision_profiles;
%             coll_local.nu_iH0 = nu_struct;
%             diff_local = compute_diffusivity_profile(velocity_results, coll_local, r_fine);
%             exp_vec = diff_local.profile_avg(:);
%             tvec = theory_primary(:);
%             valid_mask = eval_mask & isfinite(exp_vec) & isfinite(tvec) & (tvec >= 0);
%             if sum(valid_mask) < round(0.5 * sum(eval_mask))
%                 RMSE(ir0, iw, ia) = NaN;
%                 continue;
%             end
%             d = exp_vec(valid_mask) - tvec(valid_mask);
%             rmse_val = sqrt(mean(d.^2));
%             RMSE(ir0, iw, ia) = rmse_val;
%             if rmse_val < best.rmse
%                 best.rmse = rmse_val;
%                 best.amp = amp;
%                 best.width = w;
%                 best.r0 = r0;
%                 best.diff = diff_local;
%                 best.nH0 = nH0;
%                 best.nu = nu_struct;
%             end
%         end
%     end
% end
% 
% % package
% results.RMSE = RMSE;
% results.amps = amps;
% results.widths = widths;
% results.r0s = r0s;
% results.best = best;
% results.r_fine = r_fine;
% results.eval_mask = eval_mask;
% results.scaling_reference = scaling_reference;
% results.scaling_variants = scaling_variants;
% 
% % plotting (kept similar to previous)
% if opts.plot_results
%     % pick best r0 plane
%     if ~isnan(best.r0)
%         [~, ir0_best] = min(abs(r0s - best.r0));
%     else
%         ir0_best = 1;
%     end
%     figure('Name','Neutral scan RMSE (best r0)','Color','w','Position',[200 200 900 500]);
%     rmse_plane = squeeze(RMSE(ir0_best,:,:));
%     imagesc(amps, widths, rmse_plane); set(gca,'YDir','normal'); colorbar;
%     xlabel('n\_H0\_amp [m^{-3}]'); ylabel('width [m]');
%     title(sprintf('RMSE vs neutral params (r0 = %.3f m)', r0s(ir0_best)));
%     hold on;
%     plot(best.amp, best.width, 'wx', 'MarkerSize', 12, 'LineWidth', 2);
%     hold off;
% 
%     % overlay best fit
%     figure('Name','Best fit overlay','Color','w','Position',[220 220 900 520]);
%     r_norm = r_fine ./ a;
%     plot(r_norm, best.diff.profile_avg, 'k-', 'LineWidth', 2, 'DisplayName', '\chi_{eff} (best fit, exp)');
%     hold on;
%     if isfield(best.diff,'ci_lower') && isfield(best.diff,'ci_upper')
%         fill([r_norm; flipud(r_norm)], [best.diff.ci_lower; flipud(best.diff.ci_upper)], [0.85 0.85 0.85], 'FaceAlpha', 0.25, 'EdgeColor','none');
%     end
%     plot(r_norm, scaling_variants.(scaling_reference{1})(:), 'b--', 'LineWidth', 1.6, 'DisplayName', scaling_reference{1});
%     legend('Location','northwest');
%     xlabel('r/a'); ylabel('\chi [m^2/s]');
%     title('Best fit experimental vs theoretical');
%     hold off;
% end
% end


% function results = scan_neutrals_fit(velocity_results, temperature_results, ...
%     variants, constants, r_fine, amps, widths, r0s, opts)
% %SCAN_NEUTRALS_FIT Performs a grid scan of neutral Gaussian parameters.
% %   Finds the optimal parameters (amplitude, width, peak position) for the
% %   neutral density profile by minimizing the RMSE between the resulting
% %   experimental chi_phi_eff and a specified theoretical variant.
% %   The scan loops are linearized for efficient parallel execution.
% 
%     arguments
%         velocity_results (1,1) struct, temperature_results (1,1) struct
%         variants (1,1) struct, constants (1,1) struct, r_fine (:,1) double
%         amps (1,:) double, widths (1,:) double, r0s (1,:) double
%         opts (1,1) struct
%     end
% 
%     % --- 1. Setup ---
%     safe_mask = variants.meta.safe_mask;
% 
%     ref_parts = split(opts.reference_theory, '.');
%     try
%         theory_ref_struct = variants.(ref_parts{1}).(ref_parts{2});
%         theory_primary = theory_ref_struct.profile_avg;
%     catch
%         error('Invalid reference theory string: ''%s''', opts.reference_theory);
%     end
% 
%     % --- 2. Linearize Scan Parameters for Efficient Parallelization ---
%     [R0_grid, W_grid, A_grid] = ndgrid(r0s, widths, amps);
%     params_list = [R0_grid(:), W_grid(:), A_grid(:)];
%     num_calcs = size(params_list, 1);
%     rmse_vector = nan(num_calcs, 1);
% 
%     fprintf('Starting grid scan (%d total points) using parallel pool...\n', num_calcs);
% 
%     % The 'theory_primary' vector is a necessary broadcast variable, as all
%     % workers need the full profile to compute the RMSE. This is expected.
%     parfor i = 1:num_calcs
%         r0  = params_list(i, 1);
%         w   = params_list(i, 2);
%         amp = params_list(i, 3);
% 
%         % Create a temporary constants struct for this iteration
%         temp_constants = constants;
%         temp_constants.models.neutrals.n_H0_max_amp = amp;
%         temp_constants.models.neutrals.width = w;
%         temp_constants.models.neutrals.r_max_H_alpha = r0;
% 
%         % Recalculate the experimental diffusivity with these new params
%         neutral_prof = compute_neutral_density_profile(r_fine, temp_constants);
%         coll_prof = compute_collision_profiles(temperature_results, neutral_prof, r_fine);
%         diff_exp = compute_diffusivity_profile(velocity_results, coll_prof, r_fine);
% 
%         % Calculate RMSE on the safe radial window
%         residuals = diff_exp.profile_avg(safe_mask) - theory_primary(safe_mask);
%         rmse_vector(i) = sqrt(mean(residuals.^2));
%     end
% 
%     % --- 3. Find Best Result and Reshape RMSE data ---
%     RMSE_cube = reshape(rmse_vector, numel(r0s), numel(widths), numel(amps));
% 
%     [min_rmse, min_idx] = min(RMSE_cube(:));
%     [idx_r0, idx_w, idx_a] = ind2sub(size(RMSE_cube), min_idx);
% 
%     best = struct();
%     best.rmse = min_rmse;
%     best.amp = amps(idx_a);
%     best.width = widths(idx_w);
%     best.r0 = r0s(idx_r0);
% 
%     fprintf('Scan complete. Best RMSE = %.4f\n', best.rmse);
%     fprintf('  Best params: Amp=%.2e, Width=%.4f, r0=%.4f\n', best.amp, best.width, best.r0);
% 
%     % --- 4. Package Results ---
%     results.opts = opts;
%     results.RMSE_cube = RMSE_cube;
%     results.amps = amps;
%     results.widths = widths;
%     results.r0s = r0s;
%     results.best = best;
% 
%     % --- 5. Plotting ---
%     if opts.plot_results
%         fprintf('Generating study-specific plots...\n');
%         % (Plotting logic would go here)
%     end
% end