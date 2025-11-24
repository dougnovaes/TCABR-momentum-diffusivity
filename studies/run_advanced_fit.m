function run_advanced_fit(scaling_law_results, velocity_results, collision_profiles, constants, r_fine, results_dir)
%RUN_ADVANCED_FIT Performs inverse analysis to estimate neutral density.
%   This function assumes the theoretical effective diffusivity (e.g., Solomon)
%   is the "ground truth" for transport physics and optimizes a flexible
%   parametric model (Pseudo-Voigt + Offset) for the experimental diffusivity
%   to match it.
%
%   It then inverts the momentum transport equation to infer the neutral
%   particle density profile (n_H0) required to sustain this transport level.
%
%   Inputs:
%       scaling_law_results - Theoretical benchmarks (target for fitting).
%       velocity_results    - Contains the geometric factor (chi_phi_coeff).
%       collision_profiles  - Contains the CX rate coefficient profile.
%       constants, r_fine   - Standard inputs.
%       results_dir         - Directory to save plots/results.

    arguments
        scaling_law_results (1,1) struct
        velocity_results (1,1) struct
        collision_profiles (1,1) struct
        constants (1,1) struct
        r_fine (:,1) double
        results_dir (1,1) string
    end

    fprintf('\n=== Starting Advanced Inverse Analysis (Voigt Fit) ===\n');

    % --- 1. Setup Targets and Physics ---
    a = constants.machine.a;
    r_norm = r_fine / a;
    
    % Define the Target Curve (Benchmark)
    % Currently using Solomon (R0, mid-gradient) as the reference.
    % This can be changed or made an argument later.
    target_chi_profile = scaling_law_results.chi_eff_Solo; 
    
    % Define Evaluation Mask (avoid core/edge singularities if necessary)
    mask = r_norm >= 0.1 & r_norm <= 0.95;
    
    % Physics Factors for Inversion: chi_eff = [C_geo * <sv>_cx] * n_H0
    C_geo = velocity_results.chi_phi_coeff; % Scalar [m^2]
    sv_cx = collision_profiles.cx_rate.avg; % Profile [m^3/s] (ensure units!)
    
    % Check units of sv_cx. collision_profiles usually stores it in [cm^3/s] or [m^3/s]?
    % In compute_collision_profiles: all_cx_rates = ... * 1e-8 [cm^3/s]
    % And nu = rate * 1e-6 * nH0. So rate in struct is likely [cm^3/s].
    % We need everything in SI for inversion.
    % Re-deriving conversion factor just to be safe:
    % nu [s^-1] = nH0 [m^-3] * rate [m^3/s]
    % nu = nH0 * (rate_cm3s * 1e-6).
    % So the "Physics Factor" P(r) relating chi to nH0 is:
    % chi(r) = C_geo * (sv_cx_cm3s(r) * 1e-6) * nH0(r)
    % -> nH0(r) = chi(r) / (C_geo * sv_cx_cm3s(r) * 1e-6)
    
    conversion_factor = 1e-6;
    inv_physics_factor = 1 ./ (C_geo * sv_cx * conversion_factor);

    % --- 2. Configure Optimization (Pseudo-Voigt) ---
    % Model: f(x) = Offset + Amp * [ eta * L(x) + (1-eta) * G(x) ]
    % Params: [Amp, Center, Width(FWHM), Shape(eta), Offset]
    % x0 = [20, 0.85, 0.1, 0.5, 1]; % Initial guess
    
    % Lower/Upper bounds
    % Amp > 0, 0 < Center < 1, Width > 0, 0 <= eta <= 1, Offset >= 0
    lb = [0.1,  0.5,  0.01, 0,   0];
    ub = [100,  1.0,  0.5,  1,   50];
    x0 = [10,   0.85, 0.1,  0.5, 1]; 

    % Objective Function (RMSE)
    objective_fun = @(p) rmse_voigt(p, r_norm, target_chi_profile, mask);

    % Run Optimization
    options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp');
    fprintf('Optimizing Pseudo-Voigt parameters to match theoretical benchmark...\n');
    [p_opt, fval] = fmincon(objective_fun, x0, [],[],[],[], lb, ub, [], options);

    fprintf('Optimization complete. RMSE: %.4f\n', fval);
    fprintf('Parameters: Amp=%.2f, Center=%.2f, Width=%.3f, Eta=%.2f, Offset=%.2f\n', ...
        p_opt(1), p_opt(2), p_opt(3), p_opt(4), p_opt(5));

    % --- 3. Generate Optimized Profiles ---
    chi_fit = pseudo_voigt(p_opt, r_norm);
    
    % --- 4. INVERSION: Calculate Inferred Neutral Density ---
    % n_H0 = chi_fit / (C_geo * <sv> * 1e-6)
    nH0_inferred = chi_fit .* inv_physics_factor;
    
    % Sanitize inferred nH0 (remove Infs if sv_cx goes to zero, though unlikely inside mask)
    nH0_inferred(~isfinite(nH0_inferred)) = 0;
    nH0_inferred(nH0_inferred < 0) = 0;

    % --- 5. Plotting Results ---
    plot_advanced_results(r_norm, target_chi_profile, chi_fit, nH0_inferred, p_opt, results_dir);
    
    % Save results
    results_adv.p_opt = p_opt;
    results_adv.chi_fit = chi_fit;
    results_adv.nH0_inferred = nH0_inferred;
    results_adv.target = target_chi_profile;
    save(fullfile(results_dir, 'advanced_fit_results.mat'), 'results_adv');
end

% --- Helper Functions ---

function val = pseudo_voigt(p, x)
    % p = [Amp, Center, FWHM, Eta, Offset]
    A = p(1); x0 = p(2); w = p(3); eta = p(4); C = p(5);
    
    % Gaussian part
    sigma = w / (2 * sqrt(2 * log(2)));
    G = exp( - (x - x0).^2 / (2 * sigma^2) );
    
    % Lorentzian part
    gamma = w / 2;
    L = 1 ./ ( 1 + ((x - x0)/gamma).^2 );
    
    % Linear combination
    val = C + A * ( eta * L + (1 - eta) * G );
end

function err = rmse_voigt(p, x, y_target, mask)
    y_model = pseudo_voigt(p, x);
    residuals = y_model(mask) - y_target(mask);
    err = sqrt(mean(residuals.^2));
end

function plot_advanced_results(r_norm, target, fit, nH0, params, save_dir)
    % Plot 1: Fit Quality
    f1 = figure('Name', 'Advanced Fit: Diffusivity', 'WindowStyle', 'docked');
    ax1 = gca; hold on; box on; grid on;
    plot(ax1, r_norm, target, 'k--', 'LineWidth', 2, 'DisplayName', 'Theoretical Benchmark (Solomon)');
    plot(ax1, r_norm, fit, 'r-', 'LineWidth', 2.5, 'DisplayName', 'Optimised Profile (Voigt)');
    xlabel(ax1, 'r/a'); ylabel(ax1, '\chi_{\phi} [m^2/s]');
    title(ax1, 'Inverse Analysis: Fitting Effective Diffusivity');
    legend(ax1, 'Location', 'northwest');
    
    % Annotation of parameters
    txt = sprintf('Offset (Base): %.2f m^2/s\nPeak Loc: %.2f\nShape (\\eta): %.2f', params(5), params(2), params(4));
    text(ax1, 0.05, 0.5, txt, 'Units', 'normalized', 'BackgroundColor', 'w', 'EdgeColor', 'k');
    
    % Plot 2: Inferred Neutrals
    f2 = figure('Name', 'Inferred Neutral Density', 'WindowStyle', 'docked');
    ax2 = gca; hold on; box on; grid on;
    plot(ax2, r_norm, nH0, 'b-', 'LineWidth', 2.5, 'DisplayName', 'Inferred n_{H0} (from Theory)');
    xlabel(ax2, 'r/a'); ylabel(ax2, 'n_{H0} [m^{-3}]');
    title(ax2, 'Inferred Neutral Density Profile');
    legend(ax2, 'Location', 'best');
    
    % Save
    if ~exist(save_dir, 'dir'), mkdir(save_dir); end
    saveas(f1, fullfile(save_dir, 'adv_fit_diffusivity.png'));
    saveas(f2, fullfile(save_dir, 'adv_fit_neutrals.png'));
end