function plot_stage3_supporting_profiles(constants, derived_profiles, r_fine, save_dir)
%PLOT_STAGE3_SUPPORTING_PROFILES Generates consolidated plots of supporting physics profiles.
%   Creates two figures, meticulously replicating the structure from the original analysis:
%   1. Fundamental & Dimensionless Profiles: A 2x2 subplot figure showing
%      (a) Thermal Velocities, (b) Coulomb Logarithms, (c) q and s, and (d) Collisionality.
%   2. Collision Metrics: A 2x2 subplot figure replicating the layout of
%      Thesis Figure B.4, showing (a) Electron Collision Time, (b) Ion
%      Collision Time, (c) Electron Collision Frequency, and (d) Ion Collision Frequency.
%   Creates two figures, each opening as a tab in a docked window.

    arguments
        constants (1,1) struct
        derived_profiles (1,1) struct
        r_fine (:,1) double
        save_dir (1,1) string
    end

    fprintf('Plotting Stage 3: Supporting Physics Profiles...\n');
    r_norm = r_fine / constants.machine.a;
    colors = get(groot, 'DefaultAxesColorOrder');

    % --- Figure 1: Fundamental and Dimensionless Quantities ---
    fig1 = figure('Name', 'Stage 3: Fundamental Profiles', 'WindowStyle', 'docked');
    
    % Subplot (a): Thermal Velocities
    ax1a = subplot(2, 2, 1); hold(ax1a, 'on'); grid(ax1a, 'on'); box(ax1a, 'on');
    yyaxis(ax1a, 'left');
    plot(ax1a, r_norm, derived_profiles.v_th_e, '-', 'LineWidth', 3, 'Color', colors(1,:));
    ylabel(ax1a, 'Electron Thermal Velocity, $v_{th,e}$ (m/s)'); ax1a.YColor = colors(1,:); ylim(ax1a, [0, 2e7]);
    yyaxis(ax1a, 'right');
    plot(ax1a, r_norm, derived_profiles.v_th_i, '-', 'LineWidth', 3, 'Color', colors(2,:));
    ylabel(ax1a, 'Ion Thermal Velocity, $v_{th,i}$ (m/s)'); ax1a.YColor = colors(2,:); ylim(ax1a, [0, 3e5]);
    xlabel(ax1a, 'Normalised Radius ($r/a$)'); title(ax1a, '(a) Thermal Velocities'); xlim(ax1a, [0, 1]);
    set_publication_style(ax1a);
    
    % Subplot (b): Coulomb Logarithms
    ax1b = subplot(2, 2, 2); hold(ax1b, 'on'); grid(ax1b, 'on'); box(ax1b, 'on');
    plot(ax1b, r_norm, derived_profiles.ln_lambda_ee, 'LineWidth', 3, 'DisplayName', '$\lambda_{ee}$');
    plot(ax1b, r_norm, derived_profiles.ln_lambda_ei, 'LineWidth', 3, 'DisplayName', '$\lambda_{ei}$');
    plot(ax1b, r_norm, derived_profiles.ln_lambda_ii, 'LineWidth', 3, 'DisplayName', '$\lambda_{ii}$');
    xlabel(ax1b, 'Normalised Radius ($r/a$)'); ylabel(ax1b, 'Coulomb Logarithm, $\ln\Lambda$');
    title(ax1b, '(b) Coulomb Logarithms'); xlim(ax1b, [0, 1]); ylim(ax1b, [13, 17]);
    legend(ax1b, 'show', 'Location', 'best'); set_publication_style(ax1b);

    % Subplot (c): Safety Factor (q) and Magnetic Shear (s)
    ax1c = subplot(2, 2, 3); hold(ax1c, 'on'); grid(ax1c, 'on'); box(ax1c, 'on');
    
    % CORRECTION: Access flat fields directly
    plot(ax1c, r_norm, derived_profiles.q_profile, '-', 'LineWidth', 3, 'DisplayName', 'Safety Factor $q$');
    plot(ax1c, r_norm, derived_profiles.s_profile, '--', 'LineWidth', 3, 'DisplayName', 'Magnetic Shear $s$');
    
    xlabel(ax1c, 'Normalised Radius ($r/a$)'); ylabel(ax1c, 'Value'); title(ax1c, '(c) Safety Factor and Shear');
    xlim(ax1c, [0, 1]); ylim(ax1c, [0, 4]); legend(ax1c, 'show', 'Location', 'best'); set_publication_style(ax1c);

    % Subplot (d): Collisionality
    ax1d = subplot(2, 2, 4); hold(ax1d, 'on'); grid(ax1d, 'on'); box(ax1d, 'on');
    plot(ax1d, r_norm, derived_profiles.nu_star_i, '-', 'LineWidth', 3, 'DisplayName', 'Ion Collisionality ($\nu_{*i}$)');
    plot(ax1d, r_norm, derived_profiles.nu_star_e, '--', 'LineWidth', 3, 'DisplayName', 'Electron Collisionality ($\nu_{*e}$)');
    xlabel(ax1d, 'Normalised Radius ($r/a$)'); ylabel(ax1d, 'Collisionality, $\nu_*$');
    title(ax1d, '(d) Collisionality Profiles'); xlim(ax1d, [0, 1]); ylim(ax1d, [0, 4]);
    legend(ax1d, 'show', 'Location', 'best'); set_publication_style(ax1d);
    saveas(fig1, fullfile(save_dir, 'stage3_fundamental_profiles.png')); savefig(fig1, fullfile(save_dir, 'stage3_fundamental_profiles.fig'));

    % --- Figure 2: Collision Metrics (Recreates Thesis Fig B.4 Layout) ---
    fig2 = figure('Name', 'Stage 3: Collision Metrics', 'WindowStyle', 'docked');
    
    % Subplot (a): Electron Collision Time
    ax2a = subplot(2, 2, 1); hold(ax2a, 'on'); grid(ax2a, 'on'); box(ax2a, 'on');
    plot(ax2a, r_norm, derived_profiles.tau_e_calc, '--', 'LineWidth', 3.5, 'DisplayName', '$\tau_e (\lambda_{ei,calc})$');
    plot(ax2a, r_norm, derived_profiles.tau_e_15, '-', 'LineWidth', 4.5, 'DisplayName', '$\tau_e (\lambda = 15)$');
    plot(ax2a, r_norm, derived_profiles.tau_e_17, '-.', 'LineWidth', 3.5, 'DisplayName', '$\tau_e (\lambda = 17)$');
    xlabel(ax2a, '$r/a$'); ylabel(ax2a, '$\tau_e$ [s]'); title(ax2a, '(a) Electron Collision Time');
    xlim(ax2a, [0, 1]); ylim(ax2a, [0, 2e-5]); legend(ax2a, 'show', 'Location', 'best'); set_publication_style(ax2a);

    % Subplot (b): Ion Collision Time
    ax2b = subplot(2, 2, 2); hold(ax2b, 'on'); grid(ax2b, 'on'); box(ax2b, 'on');
    plot(ax2b, r_norm, derived_profiles.tau_i_calc, '--', 'LineWidth', 3.5, 'DisplayName', '$\tau_i (\lambda_{ii,calc})$');
    plot(ax2b, r_norm, derived_profiles.tau_i_15, '-', 'LineWidth', 4.5, 'DisplayName', '$\tau_i (\lambda = 15)$');
    plot(ax2b, r_norm, derived_profiles.tau_i_17, '-.', 'LineWidth', 3.5, 'DisplayName', '$\tau_i (\lambda = 17)$');
    xlabel(ax2b, '$r/a$'); ylabel(ax2b, '$\tau_i$ [s]'); title(ax2b, '(b) Ion Collision Time');
    xlim(ax2b, [0, 1]); ylim(ax2b, [0, 4e-4]); legend(ax2b, 'show', 'Location', 'best'); set_publication_style(ax2b);

    % Subplot (c): Electron Collision Frequency
    ax2c = subplot(2, 2, 3); hold(ax2c, 'on'); grid(ax2c, 'on'); box(ax2c, 'on');
    plot(ax2c, r_norm, derived_profiles.nu_ei_15 / 1e5, '-', 'LineWidth', 3.5, 'DisplayName', '$\nu_{e} (\lambda = 15)$');
    plot(ax2c, r_norm, derived_profiles.nu_ei_calc / 1e5, '--', 'LineWidth', 3, 'DisplayName', '$\nu_{e} (\lambda_{ei,calc})$');
    xlabel(ax2c, '$r/a$'); ylabel(ax2c, '$\nu_e$ [$\times10^5$ s$^{-1}$]');
    title(ax2c, '(c) Electron Collision Frequency'); xlim(ax2c, [0, 1]); ylim(ax2c, [0, 6]);
    legend(ax2c, 'show', 'Location', 'best'); set_publication_style(ax2c);

    % Subplot (d): Ion Collision Frequency
    ax2d = subplot(2, 2, 4); hold(ax2d, 'on'); grid(ax2d, 'on'); box(ax2d, 'on');
    plot(ax2d, r_norm, derived_profiles.nu_i_15 / 1e3, '-', 'LineWidth', 3.5, 'DisplayName', '$\nu_{i} (\lambda = 15)$');
    plot(ax2d, r_norm, derived_profiles.nu_i_calc / 1e3, '--', 'LineWidth', 3, 'DisplayName', '$\nu_{i} (\lambda_{ii,calc})$');
    xlabel(ax2d, '$r/a$'); ylabel(ax2d, '$\nu_i$ [$\times10^3$ s$^{-1}$]');
    title(ax2d, '(d) Ion Collision Frequency'); xlim(ax2d, [0, 1]); ylim(ax2d, [0, 16]);
    legend(ax2d, 'show', 'Location', 'best'); set_publication_style(ax2d);
    
    saveas(fig2, fullfile(save_dir, 'stage3_collision_metrics.png')); savefig(fig2, fullfile(save_dir, 'stage3_collision_metrics.fig'));
end