function plot_stage3_supporting_profiles(constants, derived_profiles, r_fine)
%PLOT_STAGE3_SUPPORTING_PROFILES Generates consolidated plots of supporting physics profiles.
%   Creates two figures:
%   1. Fundamental & Dimensionless Profiles (Thermal Vel, Coulomb Log, q, s, nu*).
%   2. Collision Metrics (Collision Times and Frequencies), replicating the
%      layout of Thesis Figure B.4.

fprintf('Plotting Stage 3: Supporting Physics Profiles...\n');

r_norm = r_fine / constants.machine.a;
default_colors = get(groot, 'DefaultAxesColorOrder'); % Get default color sequence

% --- Figure 1: Fundamental and Dimensionless Quantities ---
fig1 = figure('Name', 'Fundamental and Dimensionless Profiles', ...
              'Units', 'normalized', 'OuterPosition', [0.1 0.1 0.8 0.8]);

% Subplot (a): Thermal Velocities
ax1a = subplot(2, 2, 1);
hold(ax1a, 'on'); grid(ax1a, 'on'); box(ax1a, 'on');
yyaxis(ax1a, 'left');
plot(ax1a, r_norm, derived_profiles.v_th_e, '-', 'LineWidth', 3, 'Color', default_colors(1,:));
ylabel(ax1a, 'Electron Thermal Velocity, $v_{th,e}$ (m/s)');
ylim(ax1a, [0, 2e7]);
ax1a.YColor = default_colors(1,:);
yyaxis(ax1a, 'right');
plot(ax1a, r_norm, derived_profiles.v_th_i, '-', 'LineWidth', 3, 'Color', default_colors(2,:));
ylabel(ax1a, 'Ion Thermal Velocity, $v_{th,i}$ (m/s)');
ylim(ax1a, [0, 3e5]);
ax1a.YColor = default_colors(2,:);
xlabel(ax1a, 'Normalised Radius ($r/a$)');
title(ax1a, '(a) Thermal Velocities'); xlim(ax1a, [0, 1]);
set_publication_style(ax1a);
hold(ax1a, 'off');

% Subplot (b): Coulomb Logarithms
ax1b = subplot(2, 2, 2);
hold(ax1b, 'on'); grid(ax1b, 'on'); box(ax1b, 'on');
plot(ax1b, r_norm, derived_profiles.ln_lambda_ee, 'LineWidth', 3, 'DisplayName', '$\lambda_{ee}$', 'Color', default_colors(1,:));
plot(ax1b, r_norm, derived_profiles.ln_lambda_ei, 'LineWidth', 3, 'DisplayName', '$\lambda_{ei}$', 'Color', default_colors(2,:));
plot(ax1b, r_norm, derived_profiles.ln_lambda_ii, 'LineWidth', 3, 'DisplayName', '$\lambda_{ii}$', 'Color', default_colors(3,:));
xlabel(ax1b, 'Normalised Radius ($r/a$)'); ylabel(ax1b, 'Coulomb Logarithm, $\lambda = \ln\Lambda$');
title(ax1b, '(b) Coulomb Logarithms'); xlim(ax1b, [0, 1]); ylim(ax1b, [13, 17]); % Match thesis scale
legend(ax1b, 'show', 'Location', 'best');
set_publication_style(ax1b);
hold(ax1b, 'off');

% Subplot (c): Safety Factor (q) and Magnetic Shear (s)
ax1c = subplot(2, 2, 3);
hold(ax1c, 'on'); grid(ax1c, 'on'); box(ax1c, 'on');
plot(ax1c, r_norm, derived_profiles.q_profile, '-', 'LineWidth', 3, 'DisplayName', 'Safety Factor $q$', 'Color', default_colors(1,:));
plot(ax1c, r_norm, derived_profiles.s_profile, '--', 'LineWidth', 3, 'DisplayName', 'Magnetic Shear $s$', 'Color', default_colors(2,:));
xlabel(ax1c, 'Normalised Radius ($r/a$)'); ylabel(ax1c, '$q, s$');
title(ax1c, '(c) Safety Factor and Shear'); xlim(ax1c, [0, 1]); ylim(ax1c, [0, 4]);
legend(ax1c, 'show', 'Location', 'best');
set_publication_style(ax1c);
hold(ax1c, 'off');

% Subplot (d): Collisionality
ax1d = subplot(2, 2, 4);
hold(ax1d, 'on'); grid(ax1d, 'on'); box(ax1d, 'on');
plot(ax1d, r_norm, derived_profiles.nu_star_i, '-', 'LineWidth', 3, 'DisplayName', 'Ion Collisionality ($\nu_{*i}$)', 'Color', default_colors(1,:));
plot(ax1d, r_norm, derived_profiles.nu_star_e, '--', 'LineWidth', 3, 'DisplayName', 'Electron Collisionality ($\nu_{*e}$)', 'Color', default_colors(2,:));
xlabel(ax1d, 'Normalised Radius ($r/a$)'); ylabel(ax1d, 'Collisionality, $\nu_*$');
title(ax1d, '(d) Collisionality Profiles'); xlim(ax1d, [0, 1]); ylim(ax1d, [0, 4]);
legend(ax1d, 'show', 'Location', 'best');
set_publication_style(ax1d);
hold(ax1d, 'off');


% --- Figure 2: Collision Metrics (Recreates Thesis Fig B.4 Layout) ---
fig2 = figure('Name', 'Collision Metrics (Thesis Fig B.4 Layout)', ...
              'Units', 'normalized', 'OuterPosition', [0.1 0.1 0.8 0.8]);

% Subplot (a): Electron Collision Time
ax2a = subplot(2, 2, 1);
hold(ax2a, 'on'); grid(ax2a, 'on'); box(ax2a, 'on');
plot(ax2a, r_norm, derived_profiles.tau_e_calc, '--', 'LineWidth', 3.5, 'DisplayName', '$\tau_e (\lambda_{ei,calc})$', 'Color', default_colors(1,:));
plot(ax2a, r_norm, derived_profiles.tau_e_15, '-', 'LineWidth', 4.5, 'DisplayName', '$\tau_e (\lambda = 15)$', 'Color', default_colors(2,:));
plot(ax2a, r_norm, derived_profiles.tau_e_17, '-.', 'LineWidth', 3.5, 'DisplayName', '$\tau_e (\lambda = 17)$', 'Color', default_colors(3,:));
xlabel(ax2a, '$r/a$', 'Interpreter', 'latex');
ylabel(ax2a, '$\tau_e$ [s]', 'Interpreter', 'latex');
title(ax2a, '(a) Electron Collision Time'); xlim(ax2a, [0, 1]); ylim(ax2a, [0, 2e-5]); % Match thesis scale
legend(ax2a, 'show', 'Location', 'best');
set_publication_style(ax2a);
hold(ax2a, 'off');

% Subplot (b): Ion Collision Time
ax2b = subplot(2, 2, 2);
hold(ax2b, 'on'); grid(ax2b, 'on'); box(ax2b, 'on');
plot(ax2b, r_norm, derived_profiles.tau_i_calc, '--', 'LineWidth', 3.5, 'DisplayName', '$\tau_i (\lambda_{ii,calc})$', 'Color', default_colors(1,:));
plot(ax2b, r_norm, derived_profiles.tau_i_15, '-', 'LineWidth', 4.5, 'DisplayName', '$\tau_i (\lambda = 15)$', 'Color', default_colors(2,:));
plot(ax2b, r_norm, derived_profiles.tau_i_17, '-.', 'LineWidth', 3.5, 'DisplayName', '$\tau_i (\lambda = 17)$', 'Color', default_colors(3,:));
xlabel(ax2b, '$r/a$', 'Interpreter', 'latex');
ylabel(ax2b, '$\tau_i$ [s]', 'Interpreter', 'latex');
title(ax2b, '(b) Ion Collision Time'); xlim(ax2b, [0, 1]); ylim(ax2b, [0, 4e-4]); % Match thesis scale
legend(ax2b, 'show', 'Location', 'best');
set_publication_style(ax2b);
hold(ax2b, 'off');

% Subplot (c): Electron Collision Frequency
ax2c = subplot(2, 2, 3);
hold(ax2c, 'on'); grid(ax2c, 'on'); box(ax2c, 'on');

% Use the frequencies computed directly and those derived from tau (calc)
% Ensure fields exist, otherwise skip gracefully
if isfield(derived_profiles, 'nu_ei_15')
    plot(ax2c, r_norm, derived_profiles.nu_ei_15 / 1e5, '-', 'LineWidth', 3.5, 'DisplayName', '$\nu_{e} (\lambda = 15)$', 'Color', default_colors(2,:));
end
% plot the 'calculated' version: 1./tau_e_calc (if present)
if isfield(derived_profiles, 'tau_e_calc')
    plot(ax2c, r_norm, (1 ./ derived_profiles.tau_e_calc) / 1e5, '--', 'LineWidth', 3, 'DisplayName', '$\nu_{e} (\lambda_{ei,calc})$', 'Color', default_colors(1,:));
end

xlabel(ax2c, '$r/a$', 'Interpreter', 'latex');
ylabel(ax2c, '$\nu_e$ [$\times10^5$ s$^{-1}$]', 'Interpreter', 'latex');
title(ax2c, '(c) Electron Collision Frequency'); xlim(ax2c, [0, 1]); ylim(ax2c, [0, 6]); % adjust if required
legend(ax2c, 'show', 'Location', 'best');
set_publication_style(ax2c);
hold(ax2c, 'off');

% Subplot (d): Ion Collision Frequency
ax2d = subplot(2, 2, 4);
hold(ax2d, 'on'); grid(ax2d, 'on'); box(ax2d, 'on');

% Plot nu_i for lambda = 15 if available (use nu_i_15)
if isfield(derived_profiles, 'nu_i_15')
    plot(ax2d, r_norm, derived_profiles.nu_i_15 / 1e3, '-', 'LineWidth', 3.5, 'DisplayName', '$\nu_{i} (\lambda = 15)$', 'Color', default_colors(2,:));
end
% Plot calculated nu_i from tau_i_calc (1./tau_i_calc)
if isfield(derived_profiles, 'tau_i_calc')
    plot(ax2d, r_norm, (1 ./ derived_profiles.tau_i_calc) / 1e3, '--', 'LineWidth', 3, 'DisplayName', '$\nu_{i} (\lambda_{ii,calc})$', 'Color', default_colors(1,:));
end

xlabel(ax2d, '$r/a$', 'Interpreter', 'latex');
ylabel(ax2d, '$\nu_i$ [$\times10^3$ s$^{-1}$]', 'Interpreter', 'latex');
title(ax2d, '(d) Ion Collision Frequency'); xlim(ax2d, [0, 1]); ylim(ax2d, [0, 16]); % adjust if required
legend(ax2d, 'show', 'Location', 'best');
set_publication_style(ax2d);
hold(ax2d, 'off');


end