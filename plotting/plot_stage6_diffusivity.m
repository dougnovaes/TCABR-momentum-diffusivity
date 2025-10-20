function plot_stage6_diffusivity(effective_diffusivity_thesis, r_fine, constants)
r_norm = r_fine / constants.machine.a;
figure('Name','Stage6: Effective Diffusivity','NumberTitle','off');
hold on; grid on; box on;
if isfield(effective_diffusivity_thesis,'profile_avg')
    plot(r_norm, effective_diffusivity_thesis.profile_avg, 'b-','LineWidth',2);
end
if isfield(effective_diffusivity_thesis,'ci_lower')
    fill([r_norm; flipud(r_norm)], [effective_diffusivity_thesis.ci_lower; flipud(effective_diffusivity_thesis.ci_upper)], [0.7 0.9 0.9], 'FaceAlpha',0.25,'EdgeColor','none');
end
xlabel('r/a'); ylabel('\chi_\phi^{eff} [m^2/s]'); title('Effective diffusivity (thesis)');
xlim([0 1]);
drawnow; pause(0.05);
end
