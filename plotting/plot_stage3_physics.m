function plot_stage3_physics(~, derived_profiles, constants, r_fine)
r_norm = r_fine / constants.machine.a;
figure('Name','Stage3: Physics Profiles','NumberTitle','off');

subplot(1,2,1); hold on; grid on; box on;
if isfield(derived_profiles,'q_profile'), plot(r_norm, derived_profiles.q_profile, 'b-','LineWidth',1.8); end
if isfield(derived_profiles,'s_profile'), plot(r_norm, derived_profiles.s_profile, 'r--','LineWidth',1.6); end
xlabel('r/a'); title('q and s'); legend('q','s','Location','best'); xlim([0 1]);

subplot(1,2,2); hold on; grid on; box on;
if isfield(derived_profiles,'nu_star_i'), plot(r_norm, derived_profiles.nu_star_i, 'b-','LineWidth',1.6); end
if isfield(derived_profiles,'nu_star_e'), plot(r_norm, derived_profiles.nu_star_e, 'r--','LineWidth',1.6); end
xlabel('r/a'); title('Collisionality'); legend('\nu_* (i)','\nu_* (e)','Location','best'); xlim([0 1]);

drawnow; pause(0.05);
end
