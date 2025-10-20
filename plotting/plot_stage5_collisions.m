function plot_stage5_collisions(neutral_profile, collision_profiles, r_fine, constants)
r_norm = r_fine / constants.machine.a;

figure('Name','Stage5: Neutral & Collision','NumberTitle','off');
subplot(1,2,1); hold on; grid on; box on;
plot(r_norm, neutral_profile.n_H0, 'k-','LineWidth',1.6);
xlabel('r/a'); ylabel('n_{H0} (m^{-3})'); title('Neutral density'); xlim([0 1]);
subplot(1,2,2); hold on; grid on; box on;
if isfield(collision_profiles,'cx_rate')
    plot(r_norm, collision_profiles.cx_rate.avg, 'b-','LineWidth',1.6);
    if isfield(collision_profiles.cx_rate,'ci_lower')
        fill([r_norm; flipud(r_norm)], [collision_profiles.cx_rate.ci_lower; flipud(collision_profiles.cx_rate.ci_upper)], [0.8 0.9 1], 'FaceAlpha',0.25,'EdgeColor','none');
    end
    xlabel('r/a'); ylabel('<\sigma v>_{cx} [cm^3/s]'); title('Charge-exchange rate');
end
drawnow; pause(0.05);
end
