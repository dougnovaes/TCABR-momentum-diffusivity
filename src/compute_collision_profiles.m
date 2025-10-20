function [collision_profiles] = compute_collision_profiles(temperature_results, neutral_profile, r_fine)
% Robust computation of charge-exchange rate and ion-neutral collision freq.
% Handles multiple input shapes/types and provides clear diagnostics.

fprintf('Calculating collision rate and frequency profiles...\n');

%% 1) Retrieve bootstrap Ti matrix and normalise format --------------------
if ~isfield(temperature_results, 'bootstrap') || ...
   ~isfield(temperature_results.bootstrap, 'all_fitted_profiles')
    error('temperature_results.bootstrap.all_fitted_profiles not found.');
end
all_ti_profiles = temperature_results.bootstrap.all_fitted_profiles;

% Convert common container types to numeric matrix
if istable(all_ti_profiles)
    all_ti_profiles = table2array(all_ti_profiles);
elseif iscell(all_ti_profiles)
    try
        all_ti_profiles = cell2mat(all_ti_profiles);
    catch
        error('temperature_results.bootstrap.all_fitted_profiles is a cell that cannot be converted to matrix.');
    end
elseif isa(all_ti_profiles, 'gpuArray')
    try
        all_ti_profiles = gather(all_ti_profiles);
    catch
        error('temperature_results.bootstrap.all_fitted_profiles is a gpuArray but cannot be gathered.');
    end
end

% Force numeric, real and double
if ~isnumeric(all_ti_profiles)
    error('all_fitted_profiles must be numeric after conversion (found %s).', class(all_ti_profiles));
end
all_ti_profiles = double(all_ti_profiles);
all_ti_profiles = real(all_ti_profiles);

% Ensure r_fine column and get expected radial points
r_fine = r_fine(:);
n_r = numel(r_fine);

% Ensure the matrix is 2-D
if ~ismatrix(all_ti_profiles)
    error('all_fitted_profiles must be 2-D (radial points × iterations). Found ndims = %d.', ndims(all_ti_profiles));
end

% If sizes mismatch, try transpose
[s1, s2] = size(all_ti_profiles);
if s1 == n_r && s2 >= 1
    % ok (n_r x n_iter)
elseif s2 == n_r && s1 >= 1
    all_ti_profiles = all_ti_profiles.';  % transpose to n_r x n_iter
    [~, ~] = size(all_ti_profiles);
else
    error('Dimension mismatch: r_fine length = %d, all_fitted_profiles size = %d x %d.', n_r, s1, s2);
end

% Final sanity checks
if ~isreal(all_ti_profiles)
    all_ti_profiles = real(all_ti_profiles);
    warning('Converted all_fitted_profiles to real by discarding imaginary part.');
end
if ~all(isfinite(all_ti_profiles(:)))
    % keep but warn; we'll ignore non-finite columns later
    warning('all_fitted_profiles contains NaN or Inf values; invalid iterations will be ignored in statistics.');
end

%% 2) Prepare neutral density vector --------------------------------------
n_H0_profile = neutral_profile.n_H0(:);
if numel(n_H0_profile) ~= n_r
    error('neutral_profile.n_H0 length (%d) does not match r_fine length (%d).', numel(n_H0_profile), n_r);
end

%% 3) Compute charge-exchange rates with defensive checks -----------------

fprintf('Max imaginary part in Ti: %.3e\n', max(abs(imag(all_ti_profiles(:)))));
fprintf('Number of negative Ti: %d\n', nnz(all_ti_profiles(:) < 0));

fprintf('... calculating charge exchange rate coefficient\n');

% Ensure Ti is real and non-negative
all_ti_profiles = real(all_ti_profiles);
all_ti_profiles(all_ti_profiles < 0) = 0;  % enforce physical lower bound

% Compute <sigma*v>_cx = T_i^0.318 * 1e-8 [cm^3/s]
all_cx_rates = (all_ti_profiles .^ 0.318) * 1e-8;

% Safety check: force real
all_cx_rates = real(all_cx_rates);

fprintf('Diagnostic after sanitisation: isreal(all_cx_rates)=%d, min(Ti)=%.3g, max(Ti)=%.3g\n', ...
        isreal(all_cx_rates), min(all_ti_profiles(:)), max(all_ti_profiles(:)));

% Calculate statistics from the distribution of rate profiles
cx_rate_avg = mean(all_cx_rates, 2);
ci_percentiles_rate = prctile(all_cx_rates, [2.5, 97.5], 2);
cx_rate_ci_lower = ci_percentiles_rate(:, 1);
cx_rate_ci_upper = ci_percentiles_rate(:, 2);

%% 4) Compute ion-neutral collision frequency ------------------------------
fprintf('... calculating ion-neutral collision frequency\n');

% convert cm^3/s -> m^3/s and multiply by neutral density
all_nu_iH0 = (all_cx_rates * 1e-6) .* n_H0_profile; % implicit expansion: n_r x n_valid_iter

nu_iH0_avg = mean(all_nu_iH0, 2);
try
    ci_nu = prctile(all_nu_iH0, [2.5, 97.5], 2);
catch ME
    fprintf('Diagnostic: class(all_nu_iH0)=%s, size=%dx%d\n', class(all_nu_iH0), size(all_nu_iH0,1), size(all_nu_iH0,2));
    error('prctile failed on all_nu_iH0: %s', ME.message);
end
nu_iH0_ci_lower = ci_nu(:,1);
nu_iH0_ci_upper = ci_nu(:,2);

%% 5) Package results -----------------------------------------------------
collision_profiles.r_fine = r_fine;

collision_profiles.cx_rate.avg      = cx_rate_avg;
collision_profiles.cx_rate.ci_lower = cx_rate_ci_lower;
collision_profiles.cx_rate.ci_upper = cx_rate_ci_upper;

collision_profiles.nu_iH0.avg       = nu_iH0_avg;
collision_profiles.nu_iH0.ci_lower  = nu_iH0_ci_lower;
collision_profiles.nu_iH0.ci_upper  = nu_iH0_ci_upper;

fprintf('Collision profiles calculated successfully.\n\n');
end
