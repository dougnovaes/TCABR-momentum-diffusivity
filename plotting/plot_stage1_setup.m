function plot_stage1_setup(constants, ~, r_fine)
% Plot basic setup: base profiles and radial grid
r_norm = r_fine / constants.machine.a;
figure('Name','Stage1: Setup & Data','NumberTitle','off');
hold on; grid on; box on;
plot(r_norm, r_norm*0, 'k:'); % placeholder (optionally replace)
xlabel('Normalised radius r/a'); title('Stage 1: Data overview');
text(0.02, 0.9, sprintf('n_{points} = %d', numel(r_fine)), 'Units','normalized');
hold off;
drawnow; pause(0.05);
end
