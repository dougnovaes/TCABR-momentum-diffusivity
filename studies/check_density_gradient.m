function check_density_gradient(exp_data, constants, r_fine)
%CHECK_DENSITY_GRADIENT Compares local R/Ln vs Maslov's linear fit method.
%
%   Maslov Method: Linear fit of ln(ne) vs r over 0.2 <= r/a <= 0.8.
%   Current Method: Local gradient at r/a = 0.5.

    arguments
        exp_data (1,1) struct
        constants (1,1) struct
        r_fine (:,1) double
    end

    fprintf('\n=== Diagnosing Density Gradient (R/Ln) ===\n');

    % 1. Data Prep
    a = constants.machine.a;
    R0 = constants.machine.R0;
    r_norm = r_fine / a;
    ne = exp_data.ne_raw; % Raw density from file

    % 2. Current Method: Local Profile
    % L_ne_inv = -d(ln n)/dr
    grad_log_ne = gradient(log(ne), r_fine);
    RLn_local_profile = -R0 * grad_log_ne;
    
    % Value at r/a = 0.5
    [~, idx_05] = min(abs(r_norm - 0.5));
    RLn_at_05 = RLn_local_profile(idx_05);

    % 3. Maslov Method: Linear Fit on log(n) over 0.2 - 0.8
    mask_fit = (r_norm >= 0.2) & (r_norm <= 0.8);
    r_fit = r_fine(mask_fit);
    log_ne_fit = log(ne(mask_fit));
    
    % Polyfit degree 1: P(1)*r + P(2) -> Slope is P(1)
    % ln(n) = - (1/Ln) * r + C  => Slope = -1/Ln
    P = polyfit(r_fit, log_ne_fit, 1);
    slope = P(1); 
    Ln_inv_maslov = -slope;
    RLn_maslov = R0 * Ln_inv_maslov;

    % 4. Reporting
    fprintf('Current (Local at r/a=0.5): R/Ln = %.2f\n', RLn_at_05);
    fprintf('Maslov (Linear Fit 0.2-0.8): R/Ln = %.2f\n', RLn_maslov);
    fprintf('Difference: %.1f%%\n', 100 * abs(RLn_maslov - RLn_at_05)/RLn_maslov);

    % 5. Plotting
    f = figure('Name', 'Density Gradient Diagnosis', 'WindowStyle', 'docked');
    
    % Subplot 1: Density and Fit
    subplot(2,1,1); hold on; grid on; box on;
    plot(r_norm, log(ne), 'b-', 'LineWidth', 2, 'DisplayName', 'ln(n_e) Experimental');
    % Plot fit extension
    y_fit = polyval(P, r_fine(mask_fit));
    plot(r_norm(mask_fit), y_fit, 'r--', 'LineWidth', 2, 'DisplayName', 'Maslov Linear Fit (0.2-0.8)');
    xline(0.2, 'k:'); xline(0.8, 'k:');
    xlabel('r/a'); ylabel('ln(n_e)');
    legend('Location', 'best');
    title('Density Profile (Log Scale)');

    % Subplot 2: R/Ln Profile
    subplot(2,1,2); hold on; grid on; box on;
    plot(r_norm, RLn_local_profile, 'b-', 'LineWidth', 2, 'DisplayName', 'Local R/L_n(r)');
    yline(RLn_at_05, 'b--', 'LineWidth', 1.5, 'DisplayName', sprintf('Value at r/a=0.5 (%.2f)', RLn_at_05));
    yline(RLn_maslov, 'r-', 'LineWidth', 2, 'DisplayName', sprintf('Maslov Fit Value (%.2f)', RLn_maslov));
    
    xlabel('r/a'); ylabel('R/L_n');
    title('Normalized Density Gradient');
    legend('Location', 'best');
    ylim([0, max(RLn_local_profile)*1.2]);
    xlim([0, 1]);
end