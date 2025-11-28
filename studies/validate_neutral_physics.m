function validate_neutral_physics(constants, exp_data, derived_profiles, results_dir)
%VALIDATE_NEUTRAL_PHYSICS Performs a physics consistency check on the inferred neutral profile.
%   (Act 3 of the analysis pipeline)
%
%   1. Loads the inferred n_H0 profile from the "Advanced Fit" (Act 2).
%   2. Fits an exponential decay to the edge of the inferred profile to extract
%      the "Inferred" Penetration Length (lambda_inf).
%   3. Calculates the "Theoretical" Penetration Length (lambda_theo) based on
%      atomic physics (Ionization) using experimental Te/ne at the edge.
%   4. Compares the two to validate the physical realism of the inverse analysis.

    arguments
        constants (1,1) struct
        exp_data (1,1) struct
        derived_profiles (1,1) struct
        results_dir (1,1) string
    end

    fprintf('\n=== Starting Act 3: Physics Validation of Neutrals ===\n');

    % --- 1. Load Inferred Results (Act 2) ---
    adv_results_file = fullfile(results_dir, 'advanced_fit_results.mat');
    if ~exist(adv_results_file, 'file')
        warning('Advanced fit results not found. Run Act 2 (Advanced Inverse Analysis) first.');
        return;
    end
    load(adv_results_file, 'results_adv');
    
    % Extract profiles
    r_norm = results_adv.r_norm;
    nH0_inf = results_adv.nH0_inferred;
    a = constants.machine.a;
    r_fine = r_norm * a;

    % --- 2. Extract "Inferred" Penetration Length ---
    % We look at the edge region (e.g., 0.85 < r/a < 0.98) where the decay is exponential.
    mask_edge = (r_norm > 0.85) & (r_norm < 0.99);
    
    r_edge_fit = r_fine(mask_edge);
    n_edge_fit = nH0_inf(mask_edge);
    
    % Robust linear fit of log(n): ln(n) = A - r / lambda
    % Slope m = -1/lambda  -> lambda = -1/m
    if numel(n_edge_fit) > 5 && all(n_edge_fit > 0)
        p = polyfit(r_edge_fit, log(n_edge_fit), 1);
        lambda_inferred = -1 / p(1); % [m]
        
        % Generate the fitted exponential for plotting
        n_fit_curve = exp(polyval(p, r_edge_fit));
        
        fprintf('Inferred Penetration Length (from Slope): %.2f cm\n', lambda_inferred * 100);
    else
        warning('Could not fit exponential to edge neutrals (zeros or NaNs present).');
        lambda_inferred = NaN;
        n_fit_curve = [];
    end

    % --- 3. Calculate Theoretical Penetration Length ---
    % Formula: lambda_theo ~ v_th_n / (ne * <sv>_iz)
    % We evaluate this at a representative edge position (e.g., r/a = 0.9)
    [~, idx_eval] = min(abs(r_norm - 0.9));
    
    % Plasma parameters at evaluation point
    Te_edge = derived_profiles.Te_profile(idx_eval); % [eV]
    ne_edge = exp_data.ne_raw(idx_eval);             % [m^-3]
    
    % 3a. Neutral Thermal Velocity
    % Assumption: Neutrals are Franck-Condon atoms (~3 eV) or equilibrated with edge ions
    T_n_ev = 3.0; % [eV] Standard assumption for recycled neutrals
    m_H = constants.physics.mp;
    v_th_n = sqrt(2 * constants.physics.e_charge * T_n_ev / m_H);
    
    % 3b. Ionization Rate Coefficient <sv>_iz
    % Empirical formula for H ionization (Janev/Voronov approximation)
    % <sv>_iz [m^3/s]
    % Function of Te [eV]
    % Approximation valid for 1 < Te < 100 eV
    U = 13.6; % Ionization potential [eV]
    sv_iz = 9.7e-15 * sqrt(Te_edge) * exp(-U/Te_edge) / (1 + 0.6*sqrt(Te_edge)); % [cm^3/s]
    sv_iz = sv_iz * 1e-6; % Convert to [m^3/s]
    
    % 3c. Calculation
    freq_iz = ne_edge * sv_iz; % Ionization frequency [s^-1]
    lambda_theory = v_th_n / freq_iz; % [m]
    
    fprintf('Theoretical Penetration Length (Te=%.1f eV, ne=%.1e): %.2f cm\n', ...
        Te_edge, ne_edge, lambda_theory * 100);
    
    % --- 4. Plotting ---
    plot_physics_validation(r_norm, nH0_inf, r_edge_fit, n_fit_curve, lambda_inferred, lambda_theory, results_dir);

end

function plot_physics_validation(r_norm, nH0_inf, ~, n_fit, lambda_inf, lambda_theo, save_dir)
    % Create plot
    f = figure('Name', 'Act 3: Physics Validation', 'WindowStyle', 'docked');
    ax = gca; hold on; box on; grid on;
    
    % 1. Plot the full inferred profile
    plot(ax, r_norm, nH0_inf, 'b-', 'LineWidth', 3, 'DisplayName', 'Inferred n_{H0} (from Inverse Analysis)');
    
    % 2. Plot the exponential slope fit (if valid)
    if ~isempty(n_fit)
        % Convert r_fit (meters) back to r/a for plotting
        % Note: We assume r_norm = r_fine / a, so r_fit_norm = r_fit / a_machine
        % We need 'a' to normalize r_fit. Since we don't have 'a' in this subfunction,
        % we can infer it from the ratio of r_fit and the corresponding r_norm subset,
        % OR just plot against r_norm directly if we match indices.
        % Simpler: r_fit corresponds to mask_edge region of r_norm.
        
        % Re-finding mask to plot the segment correctly overlaid
        % (A bit hacky but effective for visualization without passing 'a' again)
        mask_plot = r_norm > 0.85 & r_norm < 0.99;
        r_plot = r_norm(mask_plot);
        % Ensure lengths match (interp if needed, but usually mask matches)
        if length(r_plot) == length(n_fit)
             plot(ax, r_plot, n_fit, 'r--', 'LineWidth', 2, ...
                 'DisplayName', sprintf('Exp. Fit (\\lambda_{inf} = %.1f cm)', lambda_inf*100));
        end
    end
    
    % Formatting
    xlabel('Normalised Radius (r/a)');
    ylabel('Neutral Density, n_{H0} [m^{-3}]');
    title('Act 3: Physical Consistency Check');
    legend('Location', 'southwest');
    xlim([0, 1]);
    
    % Add Annotation Box
    dim = [0.15 0.6 0.3 0.2];
    str = {
        '\bf{Penetration Length Check:}', ...
        sprintf('Inferred from Profile: %.1f cm', lambda_inf * 100), ...
        sprintf('Theoretical (Atomic): \\sim %.1f cm', lambda_theo * 100), ...
        ' ', ...
        '\it{(Assumes T_n \approx 3 eV)}'
    };
    annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on', 'BackgroundColor', 'w', 'EdgeColor', 'k');
    
    set_publication_style(ax);
    
    % Save
    adv_dir = fullfile(save_dir, 'advanced_studies');
    if ~exist(adv_dir, 'dir'), mkdir(adv_dir); end
    saveas(f, fullfile(adv_dir, 'act3_validation.png'));
    fprintf('Validation plot saved to %s\n', adv_dir);
end

function set_publication_style(ax)
    set(ax, 'FontSize', 16, 'LineWidth', 1.5, 'Box', 'on', 'FontName', 'Times New Roman');
    grid(ax, 'on');
end