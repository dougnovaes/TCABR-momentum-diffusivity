function plot_chi_eff_comparison(effective_thesis, variants, ~, constants)
% single-axis overlay plot with shaded experimental band and theoretical curves

a = constants.machine.a;
r_norm = variants.r_norm;
colors = get(groot,'DefaultAxesColorOrder');

figure('Name','chi_eff comparison','Color','w','Position',[200 200 900 560]);
hold on; grid on; box on;

% experimental with band
chi_exp = effective_thesis.chi_eff_profile(:);
if isfield(effective_thesis,'chi_eff_error')
    err = effective_thesis.chi_eff_error(:);
    h_exp = errorbar(r_norm, chi_exp, err, 'o', 'Color', colors(1,:), 'MarkerFaceColor', colors(1,:), 'DisplayName','Exp \chi_{eff}');
else
    % if ci available
    if isfield(effective_thesis,'ci_lower') && isfield(effective_thesis,'ci_upper')
        fill([r_norm; flipud(r_norm)], [effective_thesis.ci_lower(:); flipud(effective_thesis.ci_upper(:))], ...
            colors(1,:), 'FaceAlpha', 0.18, 'EdgeColor','none', 'DisplayName','Exp 95% CI');
    end
    h_exp = plot(r_norm, chi_exp, 'o-', 'Color', colors(1,:), 'MarkerFaceColor', colors(1,:), 'DisplayName','Exp \chi_{eff}');
end

% plot variants (choose styles)
plot(r_norm, variants.chi_eff_R0_mid, '-', 'Color', colors(2,:), 'LineWidth', 2, 'DisplayName','Theory: R_0, R/L_{V\phi}@mid');
plot(r_norm, variants.chi_eff_Rlocal_mid, '--', 'Color', colors(3,:), 'LineWidth', 1.8, 'DisplayName','Theory: R(r), R/L_{V\phi}@mid');
plot(r_norm, variants.chi_eff_R0_profile, ':', 'Color', colors(4,:), 'LineWidth', 1.8, 'DisplayName','Theory: R_0, R/L_{V\phi}(r)');
plot(r_norm, variants.chi_eff_Rlocal_profile, '-.', 'Color', colors(5,:), 'LineWidth', 1.8, 'DisplayName','Theory: R(r), R/L_{V\phi}(r)');
plot(r_norm, variants.chi_eff_peeters_mid, '-', 'Color', colors(6,:), 'LineWidth', 1.6, 'DisplayName',sprintf('Peeters-style (R/L_n@mid=%.2f)', variants.R_over_Ln_mid) );

xlabel('r/a'); ylabel('\chi_{\phi,eff} [m^2/s]','Interpreter','latex');
title('Experimental vs Theoretical \chi_{\phi,eff} (variants)','Interpreter','latex');
legend('Location','northwest'); xlim([0 1]); set(gca,'FontSize',12);

% Add small table of metrics (RMSE etc.) computed by helper
metrics = compare_metrics(chi_exp, [variants.chi_eff_R0_mid, variants.chi_eff_Rlocal_mid, variants.chi_eff_R0_profile, variants.chi_eff_Rlocal_profile]);
txt = sprintf('RMSEs: R0_mid=%.3e, Rloc_mid=%.3e, R0_prof=%.3e, Rloc_prof=%.3e', metrics.rmse);
annotation('textbox',[0.58 0.15 0.35 0.15],'String',txt,'EdgeColor','none','FontSize',10,'Interpreter','none');

hold off;

end

function metrics = compare_metrics(exp_vec, theory_matrix)
% exp_vec: n x 1, theory_matrix: n x m
m = size(theory_matrix,2);
rmse = zeros(1,m);
mean_ratio = zeros(1,m);
area_diff = zeros(1,m);
for k=1:m
    t = theory_matrix(:,k);
    valid = isfinite(exp_vec) & isfinite(t);
    e = exp_vec(valid); th = t(valid);
    rmse(k) = sqrt(mean((e - th).^2));
    mean_ratio(k) = mean(th ./ (e + eps));
    area_diff(k) = trapz(th) - trapz(e);
end
metrics.rmse = rmse;
metrics.mean_ratio = mean_ratio;
metrics.area_diff = area_diff;
end
