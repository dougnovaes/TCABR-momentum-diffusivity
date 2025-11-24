function results = scan_neutrals_fit(velocity_results, temperature_results, ...
    variants, constants, r_fine, amps, widths, r0s, opts, results_dir)
%SCAN_NEUTRALS_FIT Performs a grid scan of neutral Gaussian parameters.
%   Finds the optimal parameters (amplitude, width, peak position) for the
%   neutral density profile by minimizing the RMSE between the resulting
%   experimental chi_phi_eff and a specified theoretical variant.
%   The scan loops are linearized for efficient parallel execution.

    arguments
        velocity_results (1,1) struct
        temperature_results (1,1) struct
        variants (1,1) struct
        constants (1,1) struct
        r_fine (:,1) double
        amps (1,:) double
        widths (1,:) double
        r0s (1,:) double
        opts (1,1) struct
        results_dir (1,1) string
    end

    % --- 1. Setup ---
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

    % A tag abaixo silencia o aviso azul sobre variáveis broadcast
    parfor i = 1:num_calcs  
        r0  = params_list(i, 1); 
        w   = params_list(i, 2); %#ok<PFBNS>
        amp = params_list(i, 3);
        
        temp_constants = constants;
        temp_constants.models.neutrals.n_H0_max_amp = amp;
        temp_constants.models.neutrals.width = w;
        temp_constants.models.neutrals.r_max_H_alpha = r0;
        
        neutral_prof = compute_neutral_density_profile(r_fine, temp_constants);
        coll_prof = compute_collision_profiles(temperature_results, neutral_prof, r_fine);
        diff_exp = compute_diffusivity_profile(velocity_results, coll_prof, r_fine);
        
        residuals = diff_exp.profile_avg(safe_mask) - theory_primary(safe_mask); %#ok<PFBNS>
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
        
        try
            a = constants.machine.a;
            r_norm = r_fine / a;     % Define r_norm for plotting
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
            
            fig_overlay = figure('Name','Neutral Scan Best Fit Overlay','Color','w','Position',[220 220 900 520], 'WindowStyle', 'docked');
            ax_overlay = gca; hold(ax_overlay, 'on');
            
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