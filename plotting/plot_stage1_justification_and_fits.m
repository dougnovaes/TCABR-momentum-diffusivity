function plot_stage1_justification_and_fits(constants, exp_data, velocity_results, ...
    temperature_results, helander_results, r_fine, save_dir)
%PLOT_STAGE1_JUSTIFICATION_AND_FITS Generates plots for primary fits and physics validation.
%
%   PURPOSE:
%       Creates three key figures that establish the experimental validity of the work:
%       1. Ion Temperature Fit: Shows the canonical profile fit to experimental data.
%       2. Helander Validation: Compares the experimental toroidal velocity with the 
%          neoclassical prediction (Helander model), justifying the assumption of 
%          standard neoclassical viscosity driving the rotation against neutral friction.
%       3. Velocity Fit: Shows the polynomial and Fourier-Bessel reconstruction of the
%          velocity profile used for the spectral inversion method.
%
%   OUTPUTS:
%       Saves three .png/.fig files to the specified 'save_dir'.
%       Figures open as tabs in a single docked window for organization.

    arguments
        constants (1,1) struct
        exp_data (1,1) struct
        velocity_results (1,1) struct
        temperature_results (1,1) struct
        helander_results (1,1) struct % Updated from generic 'theoretical_models'
        r_fine (:,1) double
        save_dir (1,1) string
    end

    fprintf('Plotting Stage 1: Justification and Primary Fits...\n');
    r_norm = r_fine / constants.machine.a;

    % =========================================================================
    % FIGURE 1: ION TEMPERATURE PROFILE ANALYSIS
    % =========================================================================
    % Validates the canonical profile assumption for Ti.
    
    fig1 = figure('Name', 'Stage 1: Temp Fit', 'WindowStyle', 'docked');
    ax1 = gca; hold(ax1, 'on');
    
    % Experimental Data Points
    errorbar(ax1, exp_data.r_Ti_exp / constants.machine.a, exp_data.Ti_reconstructed_exp, ...
        exp_data.Ti_reconstructed_err_exp, 'ko', 'MarkerFaceColor', 'k', ...
        'MarkerSize', 8, 'LineWidth', 1.5, 'DisplayName', 'Experimental Data');
    
    % Canonical Fit Line
    plot(ax1, r_norm, temperature_results.profile_avg, 'r-', 'LineWidth', 4, ...
        'DisplayName', 'Canonical fit');
    
    % Confidence Interval
    fill(ax1, [r_norm; flipud(r_norm)], ...
        [temperature_results.profile_ci_lower; flipud(temperature_results.profile_ci_upper)], ...
        'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', '95\% Confidence Band');
    
    % Formatting
    xlabel(ax1, 'Normalised Radius ($r/a$)', 'Interpreter', 'latex');
    ylabel(ax1, 'Ion Temperature, $T_{i}$ (eV)', 'Interpreter', 'latex');
    title(ax1, 'Ion Temperature Profile Analysis', 'Interpreter', 'latex');
    xlim(ax1, [0, 1]); 
    legend(ax1, 'show', 'Location', 'northeast', 'Interpreter', 'latex');
    set_publication_style(ax1); 
    hold(ax1, 'off');
    
    saveas(fig1, fullfile(save_dir, 'stage1_temp_fit.png'));
    savefig(fig1, fullfile(save_dir, 'stage1_temp_fit.fig'));

    % =========================================================================
    % FIGURE 2: HELANDER MODEL JUSTIFICATION
    % =========================================================================
    % Validates that the rotation is consistent with neoclassical theory + friction.
    
    fig2 = figure('Name', 'Stage 1: Helander Model', 'WindowStyle', 'docked');
    ax2 = gca; hold(ax2, 'on');
    
    % Experimental Data (Velocity is negative in TCABR standard direction)
    % We plot -Vphi to show magnitude as positive for easier reading, or keep sign consistency.
    % Here assuming we plot magnitude or standard view.
    % Note: The input data usually comes with signs. If Vphi_exp is negative, 
    % plotting -Vphi makes it positive. Checking consistency with previous plots:
    
    h1 = errorbar(ax2, exp_data.r_Vphi_exp / constants.machine.a, -exp_data.Vphi_exp, ...
        exp_data.Vphi_err_exp, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8, ...
        'LineWidth', 1.5, 'DisplayName', 'Experimental data');
    
    % Polynomial Fit (Mean & CI)
    h2 = fill(ax2, [r_norm; flipud(r_norm)], ...
        [-velocity_results.poly_fit_ci_upper / 1000; -flipud(velocity_results.poly_fit_ci_lower / 1000)], ...
        'g', 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', '95\% Confidence interval');
    
    h3 = plot(ax2, r_norm, -velocity_results.poly_fit_avg / 1000, 'g--', 'LineWidth', 5, ...
        'DisplayName', '5th-order polynomial fit');
    
    % Theoretical Model (Helander)
    % CORRECTION: Accessing the correct field from the new compute_helander_model structure
    h4 = plot(ax2, r_norm, -helander_results.Vphi / 1000, 'r-', 'LineWidth', 4, ...
        'DisplayName', 'Helander model');
    
    yline(ax2, 0, '--', 'HandleVisibility', 'off');
    
    % Formatting
    xlabel(ax2, 'Normalised Radius ($r/a$)', 'Interpreter', 'latex'); 
    ylabel(ax2, 'Toroidal Velocity, $V_{\phi}$ (km/s)', 'Interpreter', 'latex');
    title(ax2, 'Comparison with Helander Model', 'Interpreter', 'latex');
    legend(ax2, [h1, h3, h2, h4], 'Location', 'best', 'Interpreter', 'latex');
    xlim(ax2, [0, 1]); ylim(ax2, [-10, 30]);
    set_publication_style(ax2); 
    hold(ax2, 'off');
    
    saveas(fig2, fullfile(save_dir, 'stage1_helander_justification.png'));
    savefig(fig2, fullfile(save_dir, 'stage1_helander_justification.fig'));

    % =========================================================================
    % FIGURE 3: VELOCITY PROFILE FOURIER-BESSEL ANALYSIS
    % =========================================================================
    % Shows the specific fit used for the spectral inversion method.
    
    fig3 = figure('Name', 'Stage 1: Velocity Fit', 'WindowStyle', 'docked');
    ax3 = gca; hold(ax3, 'on');
    
    errorbar(ax3, exp_data.r_Vphi_exp / constants.machine.a, -exp_data.Vphi_exp, ...
        exp_data.Vphi_err_exp, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8, ...
        'LineWidth', 1.5, 'DisplayName', 'Experimental Data');
    
    plot(ax3, r_norm, -velocity_results.poly_fit_avg / 1000, 'g--', 'LineWidth', 6, ...
        'DisplayName', '5th-order polynomial fit');
    
    fill(ax3, [r_norm; flipud(r_norm)], ...
        [-velocity_results.poly_fit_ci_lower / 1000; flipud(-velocity_results.poly_fit_ci_upper / 1000)], ...
        'g', 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', '95\% Confidence interval');
    
    plot(ax3, r_norm, -velocity_results.bessel_fit / 1000, 'r-', 'LineWidth', 4, ...
        'DisplayName', 'Fourier-Bessel series fit');
    
    yline(ax3, 0, '--', 'HandleVisibility', 'off');
    
    % Formatting
    xlabel(ax3, 'Normalised Radius ($r/a$)', 'Interpreter', 'latex'); 
    ylabel(ax3, 'Toroidal Velocity, $V_{\phi}$ (km/s)', 'Interpreter', 'latex');
    title(ax3, 'Velocity Profile Fourier-Bessel Analysis', 'Interpreter', 'latex');
    xlim(ax3, [0, 1]); ylim(ax3, [-10, 30]);
    legend(ax3, 'show', 'Location', 'best', 'Interpreter', 'latex');
    set_publication_style(ax3); 
    hold(ax3, 'off');
    
    saveas(fig3, fullfile(save_dir, 'stage1_velocity_fit.png'));
    savefig(fig3, fullfile(save_dir, 'stage1_velocity_fit.fig'));
end