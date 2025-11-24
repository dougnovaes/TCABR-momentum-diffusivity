function plot_stage2_thesis_method(constants, neutral_profile, collision_profiles, effective_diffusivity_thesis, r_fine, save_dir)
%PLOT_STAGE2_THESIS_METHOD Generates plots detailing the thesis's core method results.
%   This function creates two figures, each with two side-by-side subplots,
%   to clearly visualise the components and the final result of the thesis's
%   method for calculating effective diffusivity.
%   - Figure 1: Neutral Density (Thesis Fig 3.7 main) & CX Rate (Fig 3.7 inset).
%   - Figure 2: Collision Freq (Fig 3.8 inset) & Eff. Diffusivity (Fig 3.8 main).
%   Creates two figures, each opening as a tab in a docked window.

    arguments
        constants (1,1) struct
        neutral_profile (1,1) struct
        collision_profiles (1,1) struct
        effective_diffusivity_thesis (1,1) struct
        r_fine (:,1) double
        save_dir (1,1) string
    end

    fprintf('Plotting Stage 2: Thesis Method Results...\n');
    r_norm = r_fine / constants.machine.a;
    orange_color = [217, 83, 25] / 255;

    % --- Figure 1: Neutral Density & CX Rate (Side-by-Side) ---
    fig1 = figure('Name', 'Stage 2: Neutrals & CX Rate', 'WindowStyle', 'docked');
    
    % Subplot 1: Neutral Density Profile
    ax1_left = subplot(1, 2, 1); hold(ax1_left, 'on');
    plot(ax1_left, r_norm, neutral_profile.n_H0 / 1e16, 'LineWidth', 5, 'Color', orange_color, 'DisplayName', 'Neutral Density Profile');
    % FIX: Ensure HandleVisibility is off to hide "data1"
    yline(ax1_left, 0, '--', 'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off'); 
    grid(ax1_left, 'on');
    xlabel(ax1_left, 'Normalised Radius ($r/a$)'); ylabel(ax1_left, '$n_{H0}$ ($\times10^{16}$ m$^{-3}$)');
    title(ax1_left, 'Neutral Density Profile');
    xlim(ax1_left, [0, 1]); ylim(ax1_left, [0, 3]); 
    legend(ax1_left, 'show', 'Location', 'southeast');
    set_publication_style(ax1_left); hold(ax1_left, 'off');

    % Subplot 2: CX Rate Coefficient
    ax1_right = subplot(1, 2, 2); hold(ax1_right, 'on');
    cx_rate = collision_profiles.cx_rate;
    fill(ax1_right, [r_norm; flipud(r_norm)], [cx_rate.ci_lower; flipud(cx_rate.ci_upper)] * 1e8, ...
        orange_color, 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'DisplayName', '95\% CI');
    plot(ax1_right, r_norm, cx_rate.avg * 1e8, 'LineWidth', 3.5, 'Color', orange_color, 'DisplayName', 'Mean');
    grid(ax1_right, 'on'); xlabel(ax1_right, 'Normalised Radius ($r/a$)');
    ylabel(ax1_right, '$\langle \sigma v \rangle_{cx}$ ($\times10^{-8}$ cm$^3$/s)');
    title(ax1_right, 'CX Rate Coefficient');
    xlim(ax1_right, [0, 1]); 
    legend(ax1_right, 'show', 'Location', 'southwest'); 
    set_publication_style(ax1_right); hold(ax1_right, 'off');
    
    saveas(fig1, fullfile(save_dir, 'stage2_neutrals_cx.png')); savefig(fig1, fullfile(save_dir, 'stage2_neutrals_cx.fig'));

    % --- Figure 2: Collision Freq & Effective Diffusivity ---
    fig2 = figure('Name', 'Stage 2: Freq & Diffusivity', 'WindowStyle', 'docked');
    
    % Subplot 1: Ion-Neutral Collision Frequency
    ax2_left = subplot(1, 2, 1); hold(ax2_left, 'on');
    nu_iH0 = collision_profiles.nu_iH0;
    fill(ax2_left, [r_norm; flipud(r_norm)], [nu_iH0.ci_lower / 1e3; flipud(nu_iH0.ci_upper / 1e3)], ...
        orange_color, 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'DisplayName', '95\% CI');
    plot(ax2_left, r_norm, nu_iH0.avg / 1e3, 'LineWidth', 3.5, 'Color', orange_color, 'DisplayName', 'Mean');
    grid(ax2_left, 'on'); xlabel(ax2_left, 'Normalised Radius ($r/a$)');
    ylabel(ax2_left, 'Collision Frequency, $\nu_{iH0}$ ($\times10^3$ s$^{-1}$)');
    title(ax2_left, 'Ion-Neutral Collision Frequency');
    xlim(ax2_left, [0, 1]); 
    legend(ax2_left, 'show', 'Location', 'northwest'); 
    set_publication_style(ax2_left); hold(ax2_left, 'off');

    % Subplot 2: Effective Diffusivity
    ax2_right = subplot(1, 2, 2); hold(ax2_right, 'on');
    fill(ax2_right, [r_norm; flipud(r_norm)], [effective_diffusivity_thesis.ci_lower; flipud(effective_diffusivity_thesis.ci_upper)], ...
        [0, 0.4, 0.7], 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'DisplayName', '95\% CI');
    plot(ax2_right, r_norm, effective_diffusivity_thesis.profile_avg, 'LineWidth', 4, 'Color', [0, 0.2, 0.5], 'DisplayName', 'Mean');
    % FIX: Ensure HandleVisibility is off to hide "data1"
    yline(ax2_right, 0, '--', 'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off'); 
    grid(ax2_right, 'on');
    xlabel(ax2_right, 'Normalised Radius ($r/a$)');
    ylabel(ax2_right, 'Effective Diffusivity, $\chi_{\phi}^{\mathrm{eff}}$ [m$^2$/s]');
    title(ax2_right, 'Effective Momentum Diffusivity');
    xlim(ax2_right, [0, 1]); ylim(ax2_right, [0, 10]); 
    legend(ax2_right, 'show', 'Location', 'northwest'); 
    set_publication_style(ax2_right); hold(ax2_right, 'off');
    
    saveas(fig2, fullfile(save_dir, 'stage2_nu_diffusivity.png')); savefig(fig2, fullfile(save_dir, 'stage2_nu_diffusivity.fig'));
end