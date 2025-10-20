function plot_stage2_profiles(exp_data, velocity_results, temperature_results, r_fine, constants)
% Profile plots: velocity and ion temperature
r_norm = r_fine / constants.machine.a;

% Velocity
figure('Name','Stage2: Velocity Profile','NumberTitle','off');
hold on; grid on; box on;
errorbar(exp_data.r_Vphi_exp / constants.machine.a, -exp_data.Vphi_exp, exp_data.Vphi_err_exp, 'ko', 'MarkerFaceColor','k');
if isfield(velocity_results,'poly_fit_avg')
    plot(r_norm, -velocity_results.poly_fit_avg/1000, 'g--','LineWidth',2);
end
if isfield(velocity_results,'bessel_fit')
    plot(r_norm, -velocity_results.bessel_fit/1000, 'r-','LineWidth',1.5);
end
xlabel('r/a'); ylabel('V_\phi (km/s)'); title('Velocity profiles'); legend('Exp','Poly fit','Bessel fit','Location','best');
hold off;

% Temperature
figure('Name','Stage2: Ion Temperature','NumberTitle','off');
hold on; grid on; box on;
errorbar(exp_data.r_Ti_exp / constants.machine.a, exp_data.Ti_reconstructed_exp, exp_data.Ti_reconstructed_err_exp, 'ko', 'MarkerFaceColor','k');
if isfield(temperature_results,'profile_avg')
    plot(r_norm, temperature_results.profile_avg, 'r-','LineWidth',2);
end
if isfield(temperature_results,'profile_ci_lower') && isfield(temperature_results,'profile_ci_upper')
    fill([r_norm; flipud(r_norm)], [temperature_results.profile_ci_lower; flipud(temperature_results.profile_ci_upper)], [1 0.8 0.8], 'FaceAlpha', 0.25, 'EdgeColor','none');
end
xlabel('r/a'); ylabel('T_i (eV)'); title('Ion temperature'); legend('Reconstructed','Mean fit','95% CI','Location','best');
hold off;

drawnow; pause(0.05);
end
