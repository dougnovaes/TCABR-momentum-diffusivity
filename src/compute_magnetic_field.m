function [magnetic_field] = compute_magnetic_field(r_fine, constants)
%COMPUTE_MAGNETIC_FIELD Calculates the magnetic field profiles.
%   This function computes the toroidal, poloidal, and total magnetic
%   field profiles across the minor radius of the tokamak. The poloidal
%   field calculation is based on the model from Gilson's thesis.
%
%   Syntax:
%       magnetic_field = compute_magnetic_field(r_fine, constants)
%
%   Inputs:
%       r_fine    - The fine radial grid vector [m].
%       constants - A structure containing all machine and physical parameters.
%
%   Output:
%       magnetic_field - A structure containing the B_tor, B_pol, and B_total profiles.

fprintf('Calculating magnetic field profiles...\n');

% --- Extract necessary parameters ---
mu0 = constants.physics.mu0;
Ip  = constants.machine.Ip;
a   = constants.machine.a;
R0  = constants.machine.R0;
B0  = constants.machine.B0;

% --- Calculate Field Components ---
% Poloidal magnetic field [T] (from Gilson's thesis, Appendix A)
B_pol = (mu0 * Ip / (2 * pi * a^2)) .* r_fine .* (2 - (r_fine ./ a).^2);

% Toroidal magnetic field [T], considering toroidal effect
R_coord = R0 * (1 + r_fine / R0);
B_tor = R0 * B0 ./ R_coord;

% Total magnetic field magnitude [T]
B_total = sqrt(B_pol.^2 + B_tor.^2);

% --- Package Results into Output Structure ---
magnetic_field.r_fine   = r_fine;
magnetic_field.poloidal = B_pol;
magnetic_field.toroidal = B_tor;
magnetic_field.total    = B_total;

fprintf('Magnetic field profiles calculated successfully.\n\n');

end