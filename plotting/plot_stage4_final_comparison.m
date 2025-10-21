function plot_stage4_final_comparison(effective_diffusivity_thesis, scaling_law_results, constants, r_fine)
% Simple plot: overlays experimental effective diffusivity and theory curves.
% Assumes fields exist:
% - effective_diffusivity_thesis.profile_avg (or .chi_eff_profile)
% - scaling_law_results.chi_phi_solo
% - scaling_law_results.chi_eff_Solo, _Hahm, _Gurcan, _Peeters_Rln2, _Peeters_Rln_calc

% Normalise radius
a = constants.machine.a;
r_norm = r_fine(:) / a;

% Fetch experimental profile (support both names)
if isfield(effective_diffusivity_thesis,'chi_eff_profile')
    chi_exp = effective_diffusivity_thesis.chi_eff_profile(:);
elseif isfield(effective_diffusivity_thesis,'profile_avg')
    chi_exp = effective_diffusivity_thesis.profile_avg(:);
else
    chi_exp = [];
end

% Theoretical curves (assume existence)
chi_phi_solo = scaling_law_results.chi_phi_solo(:);
chi_eff_Solo  = scaling_law_results.chi_eff_Solo(:);
chi_eff_Hahm  = scaling_law_results.chi_eff_Hahm(:);
chi_eff_Gurcan= scaling_law_results.chi_eff_Gurcan(:);
chi_eff_Peet2 = scaling_law_results.chi_eff_Peeters_Rln2(:);
chi_eff_PeetC = scaling_law_results.chi_eff_Peeters_Rln_calc(:);

% Basic plotting (single axes)
figure('Name','Final Comparison: Experimental vs Theoretical','NumberTitle','off','Color','w','Position',[100 100 900 600]);
hold on; grid on; box on;

% Colours (default cycle)
cols = get(groot,'DefaultAxesColorOrder');

% Plot experimental (if present)
if ~isempty(chi_exp)
    plot(r_norm, chi_exp, 'o-','Color',cols(1,:),'MarkerFaceColor',cols(1,:),'LineWidth',1.4,'DisplayName','\chi_{eff}^{(Exp)}');
end

% Plot theory
plot(r_norm, chi_phi_solo, '-', 'Color', cols(2,:), 'LineWidth', 1.8, 'DisplayName','\chi_{\phi}^{(Solo)}');
plot(r_norm, chi_eff_Solo,  '--','Color', cols(3,:), 'LineWidth', 1.6, 'DisplayName','\chi_{eff}^{(Solo)}');
plot(r_norm, chi_eff_Hahm,  '-.','Color', cols(4,:), 'LineWidth', 1.6, 'DisplayName','\chi_{eff}^{(Hahm)}');
plot(r_norm, chi_eff_Gurcan,':','Color', cols(5,:), 'LineWidth', 1.6, 'DisplayName','\chi_{eff}^{(Gurcan)}');
plot(r_norm, chi_eff_Peet2, '-', 'Color', cols(6,:), 'LineWidth', 1.4, 'DisplayName','\chi_{eff}^{(Peeters, R/L_n=2)}');
plot(r_norm, chi_eff_PeetC, '--','Color', cols(7,:), 'LineWidth', 1.4, 'DisplayName',sprintf('\\chi_{eff}^{(Peeters, R/L_n=%.2f)}', scaling_law_results.R_over_Ln_const));

xlabel('r/a','Interpreter','latex'); 
ylabel('\chi_{\phi} [m^2/s]','Interpreter','latex');
title('Experimental vs Theoretical Effective Momentum Diffusivity','Interpreter','latex');
legend('Location','northwest','Interpreter','latex');
xlim([0 1]); 
set(gca,'FontSize',12);

% Save figure
save_dir = fullfile('results','plots');
if ~exist(save_dir,'dir'), mkdir(save_dir); end
saveas(gcf, fullfile(save_dir,'stage4_final_comparison_simple.png'));

end
