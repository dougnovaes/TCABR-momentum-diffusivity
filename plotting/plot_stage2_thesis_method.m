function plot_stage2_thesis_method(constants, neutral_profile, collision_profiles, effective_diffusivity_thesis, r_fine)
%PLOT_STAGE2_THESIS_METHOD Generates plots detailing the thesis's core method results.
%   This function creates two figures, each with two side-by-side subplots,
%   to clearly visualise the components and the final result of the thesis's
%   method for calculating effective diffusivity. This layout replaces the
%   previous inset plot style as requested.

    arguments
        constants (1,1) struct, neutral_profile (1,1) struct, collision_profiles (1,1) struct
        effective_diffusivity_thesis (1,1) struct, r_fine (:,1) double, save_dir (1,1) string
    end

    fprintf('Plotting Stage 2: Thesis Method Results (Side-by-Side Layout)...\n');
    
    r_norm = r_fine / constants.machine.a;
    orange_color = [217, 83, 25] / 255; % '#D95319'
    
    % --- Figure 1: Neutral Density & CX Rate (Side-by-Side) ---
    fig1 = figure('Name', 'Neutral Density and CX Rate', 'Position', ...
        [100, 100, 1200, 600]); % Adjust size/position
    fig1.Position(3) = fig1.Position(3) * 1; % Widen the figure
    
    % Subplot 1: Neutral Density Profile (Matches Thesis Fig. 3.7 Main)
    ax1_left = subplot(1, 2, 1);
    hold(ax1_left, 'on');
    plot(ax1_left, r_norm, neutral_profile.n_H0 / 1e16, 'LineWidth', 5, 'Color', orange_color, ...
         'DisplayName', 'Neutral particles density');
    yline(ax1_left, 0, '--', 'HandleVisibility', 'off');
    grid(ax1_left, 'on');
    xlabel(ax1_left, 'Normalised Radius ($r/a$)');
    ylabel(ax1_left, '$n_{H0}$ ($\times10^{16}$ m$^{-3}$)', 'Interpreter', 'latex');
    title(ax1_left, 'Neutral Density Profile');
    xlim(ax1_left, [0, 1]); ylim(ax1_left, [0, 3]);
    xticks(ax1_left, 0:0.2:1); yticks(ax1_left, 0:1:3);
    legend(ax1_left, 'show', 'Location', 'southeast');
    set_publication_style(ax1_left);
    hold(ax1_left, 'off');
    
    % Subplot 2: CX Rate Coefficient (Matches Thesis Fig. 3.7 Inset)
    ax1_right = subplot(1, 2, 2);
    hold(ax1_right, 'on');
    cx_rate = collision_profiles.cx_rate;
    fill(ax1_right, [r_norm; flipud(r_norm)], [cx_rate.ci_lower; flipud(cx_rate.ci_upper)] * 1e8, ...
        orange_color, 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'DisplayName', '95% Confidence band');
    plot(ax1_right, r_norm, cx_rate.avg * 1e8, 'LineWidth', 3.5, 'Color', orange_color, ...
         'DisplayName', 'Reaction rate coefficient');
    grid(ax1_right, 'on');
    xlabel(ax1_right, 'Normalised Radius ($r/a$)');
    ylabel(ax1_right, '$\langle \sigma v \rangle_{cx}$ ($\times10^{-8}$ cm$^3$/s)', 'Interpreter', 'latex');
    title(ax1_right, 'CX Rate Coefficient');
    xlim(ax1_right, [0, 1]); ylim(ax1_right, [2, 7]);
    xticks(ax1_right, 0:0.2:1); yticks(ax1_right, 2:1:7);
    legend(ax1_right, 'show', 'Location', 'southwest');
    set_publication_style(ax1_right);
    hold(ax1_right, 'off');
    
    % --- Figure 2: Effective Diffusivity & Collision Freq. (Side-by-Side) ---
    fig2 = figure('Name', 'Effective Diffusivity and Collision Frequency (Side-by-Side)', ...
                  'Units', 'normalized', 'OuterPosition', [0.1 0.1 .8 .8]);
    fig2.Position(3) = fig2.Position(3) * 1; % Widen the figure
    
    % Subplot 1: Ion-Neutral Collision Frequency (Matches Thesis Fig. 3.8 Inset)
    ax2_right = subplot(1, 2, 1);
    hold(ax2_right, 'on');
    nu_iH0 = collision_profiles.nu_iH0;
    fill(ax2_right, [r_norm; flipud(r_norm)], [nu_iH0.ci_lower / 1e3; flipud(nu_iH0.ci_upper / 1e3)], ...
        orange_color, 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'DisplayName', '95% Confidence band');
    plot(ax2_right, r_norm, nu_iH0.avg / 1e3, 'LineWidth', 3.5, 'Color', orange_color, ...
         'DisplayName', 'Collision frequency');
    grid(ax2_right, 'on');
    xlabel(ax2_right, 'Normalised Radius ($r/a$)');
    ylabel(ax2_right, 'Collision Frequency, $\nu_{iH0}$ ($\times10^3$ s$^{-1}$)');
    title(ax2_right, 'Ion-Neutral Collision Frequency');
    xlim(ax2_right, [0, 1]); ylim(ax2_right, [0, 1]);
    xticks(ax2_right, 0:0.2:1); yticks(ax2_right, 0:0.5:1);
    legend(ax2_right, 'show', 'Location', 'northwest');
    set_publication_style(ax2_right);
    hold(ax2_right, 'off');
    
    % Subplot 2: Effective Diffusivity (Matches Thesis Fig. 3.8 Main)
    ax2_left = subplot(1, 2, 2);
    hold(ax2_left, 'on');
    eff_thesis = effective_diffusivity_thesis;
    fill(ax2_left, [r_norm; flipud(r_norm)], [eff_thesis.ci_lower; flipud(eff_thesis.ci_upper)], ...
        [0, 222/255, 230/255], 'FaceAlpha', 0.55, 'EdgeColor', 'none', 'DisplayName', '95% Confidence Band');
    plot(ax2_left, r_norm, eff_thesis.profile_avg, 'LineWidth', 4, 'Color', [0, 0.2470, 0.5410], ...
         'DisplayName', 'Mean effective diffusivity');
    yline(ax2_left, 0, '--', 'HandleVisibility', 'off');
    grid(ax2_left, 'on');
    xlabel(ax2_left, 'Normalised Radius ($r/a$)');
    ylabel(ax2_left, 'Effective Diffusivity, $\chi_{\phi}^{eff}$ [m$^2$/s]');
    title(ax2_left, 'Effective Momentum Diffusivity');
    xlim(ax2_left, [0, 1]); ylim(ax2_left, [0, 10]);
    xticks(ax2_left, 0:0.2:1); yticks(ax2_left, 0:5:10);
    legend(ax2_left, 'show', 'Location', 'northwest');
    set_publication_style(ax2_left);
    hold(ax2_left, 'off');

end