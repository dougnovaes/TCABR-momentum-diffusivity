function plot_stage4_models(theoretical_models, derived_profiles, constants, r_fine)
r_norm = r_fine / constants.machine.a;

% Helander compare
figure('Name','Stage4: Helander vs Fit','NumberTitle','off');
hold on; grid on; box on;
if isfield(theoretical_models,'Vphi_Helander')
    plot(r_norm, -theoretical_models.Vphi_Helander/1000, 'r-', 'LineWidth', 2, 'DisplayName', 'Helander');
end
if isfield(derived_profiles,'Vphi_exp_fit') % optional
    plot(r_norm, -derived_profiles.Vphi_exp_fit/1000, 'g--', 'LineWidth',1.5, 'DisplayName','Fit');
end
xlabel('r/a'); ylabel('V_\phi (km/s)'); title('Helander model');
legend('Location','best');
hold off;

% Solomon diffusivity plot (if present)
if isfield(theoretical_models,'chi_phi_solo_Rln_calc')
    figure('Name','Stage4: Solomon Diffusivity','NumberTitle','off');
    chi = theoretical_models.chi_phi_solo_Rln_calc;
    plot(r_norm, chi.profile, 'r-','LineWidth',1.8); hold on;
    if isfield(chi,'lower_bound') && isfield(chi,'upper_bound')
        fill([r_norm; flipud(r_norm)], [chi.lower_bound; flipud(chi.upper_bound)], [0.85 0.325 0.098], 'FaceAlpha', 0.2, 'EdgeColor','none');
    end
    xlabel('r/a'); ylabel('\chi_\phi (m^2/s)'); title('Solomon model (R/L_n ~ mid)');
    hold off;
end

drawnow; pause(0.05);
end
