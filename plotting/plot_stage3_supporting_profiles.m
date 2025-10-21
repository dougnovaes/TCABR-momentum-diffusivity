function plot_stage3_supporting_profiles(constants, derived_profiles, r_fine)
%PLOT_STAGE3_SUPPORTING_PROFILES Generates consolidated plots of supporting physics profiles.
%   This function creates two figures, each containing multiple subplots,
%   to visualise key derived profiles such as thermal velocities, collision
%   metrics, safety factor, shear, and collisionality, matching the layout
%   agreed upon in the plotting plan.

fprintf('Plotting Stage 3: Supporting Physics Profiles...\n');

r_norm = r_fine / constants.machine.a;

% --- Figure 1: Fundamental and Dimensionless Quantities ---
fig1 = figure('Name', 'Fundamental and Dimensionless Profiles', ...
              'Units', 'normalized', 'OuterPosition', [0 0.05 0.8 0.9]);

% Subplot (a): Thermal Velocities
ax1a = subplot(2, 2, 1);
hold(ax1a, 'on'); grid(ax1a, 'on'); box(ax1a, 'on');
% NOTE: 'colororder' line removed. yyaxis handles color cycling.
yyaxis(ax1a, 'left');
plot(ax1a, r_norm, derived_profiles.v_th_e, '-', 'LineWidth', 3);
ylabel(ax1a, 'Electron Thermal Velocity, $v_{th,e}$ (m/s)');
ylim(ax1a, [0, 2e7]);
ax1a.YColor = [0 0.4470 0.7410]; % Blue
yyaxis(ax1a, 'right');
plot(ax1a, r_norm, derived_profiles.v_th_i, '-', 'LineWidth', 3);
ylabel(ax1a, 'Ion Thermal Velocity, $v_{th,i}$ (m/s)');
ylim(ax1a, [0, 3e5]);
ax1a.YColor = [0.8500 0.3250 0.0980]; % Orange
xlabel(ax1a, 'Normalised Radius ($r/a$)');
title(ax1a, '(a) Thermal Velocities'); xlim(ax1a, [0, 1]);
set_publication_style(ax1a);
hold(ax1a, 'off');

% Subplot (b): Coulomb Logarithms
ax1b = subplot(2, 2, 2);
hold(ax1b, 'on'); grid(ax1b, 'on'); box(ax1b, 'on');
plot(ax1b, r_norm, derived_profiles.ln_lambda_ee, 'LineWidth', 3, 'DisplayName', '$\lambda_{ee}$');
plot(ax1b, r_norm, derived_profiles.ln_lambda_ei, 'LineWidth', 3, 'DisplayName', '$\lambda_{ei}$');
plot(ax1b, r_norm, derived_profiles.ln_lambda_ii, 'LineWidth', 3, 'DisplayName', '$\lambda_{ii}$');
xlabel(ax1b, 'Normalised Radius ($r/a$)'); ylabel(ax1b, 'Coulomb Logarithm, $\lambda = \ln\Lambda$');
title(ax1b, '(b) Coulomb Logarithms'); xlim(ax1b, [0, 1]);
legend(ax1b, 'show', 'Location', 'best');
set_publication_style(ax1b);
hold(ax1b, 'off');

% Subplot (c): Safety Factor (q) and Magnetic Shear (s)
ax1c = subplot(2, 2, 3);
hold(ax1c, 'on'); grid(ax1c, 'on'); box(ax1c, 'on');
plot(ax1c, r_norm, derived_profiles.q_profile, 'b-', 'LineWidth', 3, 'DisplayName', 'Safety Factor $q$');
plot(ax1c, r_norm, derived_profiles.s_profile, 'r--', 'LineWidth', 3, 'DisplayName', 'Magnetic Shear $s$');
xlabel(ax1c, 'Normalised Radius ($r/a$)'); ylabel(ax1c, 'Value');
title(ax1c, '(c) Safety Factor and Shear'); xlim(ax1c, [0, 1]); ylim(ax1c, [0, 4]);
legend(ax1c, 'show', 'Location', 'best');
set_publication_style(ax1c);
hold(ax1c, 'off');

% Subplot (d): Collisionality
ax1d = subplot(2, 2, 4);
hold(ax1d, 'on'); grid(ax1d, 'on'); box(ax1d, 'on');
plot(ax1d, r_norm, derived_profiles.nu_star_i, 'b-', 'LineWidth', 3, 'DisplayName', 'Ion Collisionality ($\nu_{*i}$)');
plot(ax1d, r_norm, derived_profiles.nu_star_e, 'r--', 'LineWidth', 3, 'DisplayName', 'Electron Collisionality ($\nu_{*e}$)');
xlabel(ax1d, 'Normalised Radius ($r/a$)'); ylabel(ax1d, 'Collisionality, $\nu_*$');
title(ax1d, '(d) Collisionality Profiles'); xlim(ax1d, [0, 1]); ylim(ax1d, [0, 4]);
legend(ax1d, 'show', 'Location', 'best');
set_publication_style(ax1d);
hold(ax1d, 'off');


% --- Figure 2: Collision Metrics Compared ---
fig2 = figure('Name', 'Collision Metrics Comparison', ...
              'Units', 'normalized', 'OuterPosition', [0.1 0.1 0.7 0.8]);
fig2.Position(3) = fig2.Position(3) * 1.8; % Widen for side-by-side plots

% Subplot (a): Collision Times
ax2a = subplot(1, 2, 1);
hold(ax2a, 'on'); grid(ax2a, 'on'); box(ax2a, 'on');
% NOTE: 'colororder' line removed. yyaxis handles color cycling.
yyaxis(ax2a, 'left');
plot(ax2a, r_norm, derived_profiles.tau_e_calc, '--', 'LineWidth', 2.5, 'DisplayName', '$\tau_e (\lambda_{ei,calc})$');
plot(ax2a, r_norm, derived_profiles.tau_e_15, '-', 'LineWidth', 3, 'DisplayName', '$\tau_e (\lambda = 15)$');
plot(ax2a, r_norm, derived_profiles.tau_e_17, '-.', 'LineWidth', 2.5, 'DisplayName', '$\tau_e (\lambda = 17)$');
ylabel(ax2a, 'Electron Collision Time, $\tau_e$ [s]');
ax2a.YColor = [0 0.4470 0.7410]; % Blue tones
yyaxis(ax2a, 'right');
plot(ax2a, r_norm, derived_profiles.tau_i_calc, '--', 'LineWidth', 2.5, 'DisplayName', '$\tau_i (\lambda_{ii,calc})$');
plot(ax2a, r_norm, derived_profiles.tau_i_15, '-', 'LineWidth', 3, 'DisplayName', '$\tau_i (\lambda = 15)$');
plot(ax2a, r_norm, derived_profiles.tau_i_17, '-.', 'LineWidth', 2.5, 'DisplayName', '$\tau_i (\lambda = 17)$');
ylabel(ax2a, 'Ion Collision Time, $\tau_i$ [s]');
ax2a.YColor = [0.8500 0.3250 0.0980]; % Orange tones
xlabel(ax2a, 'Normalised Radius ($r/a$)');
title(ax2a, '(a) Collision Times'); xlim(ax2a, [0, 1]);
legend(ax2a, 'show', 'Location', 'best');
set_publication_style(ax2a);
hold(ax2a, 'off');

% Subplot (b): Collision Frequencies
ax2b = subplot(1, 2, 2);
hold(ax2b, 'on'); grid(ax2b, 'on'); box(ax2b, 'on');
% NOTE: 'colororder' line removed. yyaxis handles color cycling.
yyaxis(ax2b, 'left');
plot(ax2b, r_norm, derived_profiles.nu_ei / 1e5, '-', 'LineWidth', 3);
ylabel(ax2b, 'Electron Collision Frequency, $\nu_{ei}$ [$\times10^5$ s$^{-1}$]');
ax2b.YColor = [0 0.4470 0.7410]; % Blue
yyaxis(ax2b, 'right');
plot(ax2b, r_norm, derived_profiles.nu_i / 1e3, '-', 'LineWidth', 3);
ylabel(ax2b, 'Ion Collision Frequency, $\nu_i$ [$\times10^3$ s$^{-1}$]');
ax2b.YColor = [0.8500 0.3250 0.0980]; % Orange
xlabel(ax2b, 'Normalised Radius ($r/a$)');
title(ax2b, '(b) Collision Frequencies'); xlim(ax2b, [0, 1]);
set_publication_style(ax2b);
hold(ax2b, 'off');

end