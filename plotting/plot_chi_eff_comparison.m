function plot_chi_eff_comparison(diffusivity_results, variants, ~, constants, varargin)
% PLOT_CHI_EFF_COMPARISON Plot experimental chi_eff and multiple theoretical variants.
%   plot_chi_eff_comparison(diffusivity_results, variants, scaling_law_results, constants)
%
% Inputs:
%   diffusivity_results - struct with fields .profile_avg, .ci_lower, .ci_upper, .r_fine
%   variants            - output of compute_chi_eff_variants
%   scaling_law_results - (optional) original scaling results for extra info
%   constants           - to extract a, R0, etc.

r = diffusivity_results.r_fine(:);
r_norm = r ./ constants.machine.a;
colors = get(groot,'DefaultAxesColorOrder');

figure('Name','χ_φ,eff: Experimental vs Theoretical','Color','w','Position',[100 100 900 600]);
hold on; grid on; box on;

% Experimental
plot(r_norm, diffusivity_results.profile_avg, '-o', 'LineWidth', 2, 'MarkerSize', 4, 'Color', colors(1,:));
if isfield(diffusivity_results,'ci_lower') && isfield(diffusivity_results,'ci_upper')
    fill([r_norm; flipud(r_norm)], [diffusivity_results.ci_lower; flipud(diffusivity_results.ci_upper)], ...
        colors(1,:), 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility','off');
end

% Plot chi_phi_solo baseline
plot(r_norm, variants.chi_phi_solo, '-', 'LineWidth', 1.8,'Color', colors(2,:), 'DisplayName','\chi_\phi^{(Solo)} (baseline)');

% Plot R0_mid variants
fields_mid = fieldnames(variants.chi_eff_R0_mid);
for k=1:numel(fields_mid)
    fld = fields_mid{k};
    data = variants.chi_eff_R0_mid.(fld);
    plot(r_norm, data, '--', 'LineWidth', 1.6, 'DisplayName', sprintf('\\chi_{eff}^{%s} (R0, mid)', fld));
end

% Plot R0_profile variants (dashed-dotted)
fields_prof = fieldnames(variants.chi_eff_R0_profile);
for k=1:numel(fields_prof)
    fld = fields_prof{k};
    data = variants.chi_eff_R0_profile.(fld);
    plot(r_norm, data, '-.', 'LineWidth', 1.4, 'DisplayName', sprintf('\\chi_{eff}^{%s} (R0, profile)', fld));
end

% R_local variants
if isfield(variants,'chi_eff_Rlocal_profile')
    fields_rl = fieldnames(variants.chi_eff_Rlocal_profile);
    for k=1:numel(fields_rl)
        fld = fields_rl{k};
        data = variants.chi_eff_Rlocal_profile.(fld);
        plot(r_norm, data, ':', 'LineWidth', 1.6, 'DisplayName', sprintf('\\chi_{eff}^{%s} (R_{local}, profile)', fld));
    end
end

% Labels & legend
xlabel('r/a'); ylabel('\chi_{\phi,eff} [m^2/s]');
title('Experimental \chi_{\phi,eff} and theoretical variants','Interpreter','latex');
legend('Location','northwest','Interpreter','tex');
xlim([0,1]); ylim('auto');
% Annotate important scalars
xtext = 0.02; ytext = 0.95 * ylim;
if isfield(variants,'meta')
    txt = sprintf('R0=%.3f m, R/L_n(mid)=%.2f, R/L_Vphi(mid)=%.2f', variants.meta.R0, variants.meta.R_over_Ln_mid, variants.meta.R_over_LVphi_mid);
    text(xtext, ytext(2), txt, 'Units','normalized', 'FontSize', 9, 'BackgroundColor','w', 'EdgeColor','k');
end

% Save
savepath = fullfile('results','plots');
if ~exist(savepath,'dir'), mkdir(savepath); end
saveas(gcf, fullfile(savepath,'chi_eff_comparison.png'));
hold off;
end
