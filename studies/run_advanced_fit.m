function run_advanced_fit(velocity_results, collision_profiles, variants, constants, r_fine, results_dir)
%RUN_ADVANCED_FIT Performs inverse analysis to estimate neutral density.
%   Optimizes a Pseudo-Voigt model for effective diffusivity to match
%   the theoretical benchmark.

    arguments
        velocity_results (1,1) struct
        collision_profiles (1,1) struct
        variants (1,1) struct
        constants (1,1) struct
        r_fine (:,1) double
        results_dir (1,1) string
    end

    fprintf('\n=== Starting Advanced Inverse Analysis (Voigt Fit) ===\n');

    % --- 1. Setup Targets and Physics ---
    a = constants.machine.a;
    r_norm = r_fine / a;
    
    % Benchmark: Solomon (R0) using LOCAL Profile Gradients
    if isfield(scaling_law_results, 'chi_eff') && isfield(scaling_law_results.chi_eff, 'Solo')
        target_chi_profile = scaling_law_results.chi_eff.Solo;
    elseif isfield(scaling_law_results, 'chi_eff_Solo')
        target_chi_profile = scaling_law_results.chi_eff_Solo;
    else
        error('Could not find Solomon profile in scaling_law_results (checked .chi_eff.Solo and .chi_eff_Solo).');
    end
    
    target_chi_profile(isnan(target_chi_profile)) = 0;
    target_chi_profile(target_chi_profile < 0) = 0;
    
    % --- 2. Physical Validity Mask ---
    mask_range = (r_norm > 0.20) & (r_norm < 0.98);
    
    if isfield(variants, 'meta') && isfield(variants.meta, 'R_over_LVphi_profile')
        R_over_LVphi = variants.meta.R_over_LVphi_profile;
        mask_physics = abs(R_over_LVphi) > 1.0; % Filter singularities
    else
        mask_physics = true(size(r_norm));
    end
    mask = mask_range & mask_physics;
    
    fprintf('Mask applied: %d of %d points used (%.1f%%).\n', ...
        sum(mask), numel(mask), 100*sum(mask)/numel(mask));

    % Physics Factors
    C_geo = velocity_results.chi_phi_coeff; 
    sv_cx = collision_profiles.cx_rate.avg; 
    conversion_factor = 1e-6; 
    physics_factor = C_geo .* (sv_cx .* conversion_factor);
    physics_factor(physics_factor == 0) = eps;
    inv_physics_factor = 1 ./ physics_factor;

    % --- 3. Optimization (Pseudo-Voigt) ---
    % Params: [Amp, Center, Width, Eta, Offset]
    lb = [0.1,  0.5,  0.02, 0,   0];
    ub = [100,  1.0,  0.4,  1,   20];
    x0 = [10,   0.85, 0.15, 0.5, 1]; 

    objective_fun = @(p) rmse_voigt(p, r_norm, target_chi_profile, mask);
    options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp');
    
    try
        fprintf('Optimizing Pseudo-Voigt parameters...\n');
        [p_opt, fval] = fmincon(objective_fun, x0, [],[],[],[], lb, ub, [], options);
    catch ME
        warning(ME.identifier, 'Optimization failed: %s', ME.message);
        return;
    end

    fprintf('Optimization complete. RMSE: %.4f\n', fval);

    % --- 4. Inversion ---
    chi_fit = pseudo_voigt(p_opt, r_norm);
    nH0_inferred = chi_fit .* inv_physics_factor;
    nH0_inferred(~isfinite(nH0_inferred)) = 0;
    nH0_inferred(nH0_inferred < 0) = 0;

    % --- 5. Save Results ---
    results_adv.p_opt = p_opt;
    results_adv.chi_fit = chi_fit;
    results_adv.nH0_inferred = nH0_inferred;
    results_adv.target = target_chi_profile;
    results_adv.mask = mask;
    results_adv.r_norm = r_norm;
    
    % Ensure directory exists before saving (Fixing the error)
    if ~exist(results_dir, 'dir'), mkdir(results_dir); end
    
    save_path = fullfile(results_dir, 'advanced_fit_results.mat');
    save(save_path, 'results_adv');
    fprintf('Results saved to: %s\n', save_path);
    
    % --- 6. Plotting (Delegated) ---
    if exist('plot_advanced_fit_results', 'file') == 2
        plot_advanced_fit_results(results_adv, results_dir);
    else
        warning('Plotting function not found.');
    end
end

% --- Helpers ---
function val = pseudo_voigt(p, x)
    A = p(1); x0 = p(2); w = p(3); eta = p(4); C = p(5);
    sigma = w / (2 * sqrt(2 * log(2)));
    G = exp( - (x - x0).^2 / (2 * sigma^2) );
    gamma = w / 2;
    L = 1 ./ ( 1 + ((x - x0)/gamma).^2 );
    val = C + A * ( eta * L + (1 - eta) * G );
end

function err = rmse_voigt(p, x, y_target, mask)
    y_model = pseudo_voigt(p, x);
    residuals = y_model(mask) - y_target(mask);
    err = sqrt(mean(residuals.^2));
end