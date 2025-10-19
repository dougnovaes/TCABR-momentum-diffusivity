function zeros_bessel = findBesselZeros(n, order)
%FINDBESSELZEROS Finds the first n zeros of a Bessel function of a given order.
%   This function is a utility to calculate the zeros required for the
%   Fourier-Bessel series expansion.
%
%   Syntax:
%       zeros_bessel = findBesselZeros(n, order)
%
%   Inputs:
%       n     - The number of zeros to find.
%       order - The order of the Bessel function (e.g., 0 for J0).
%
%   Output:
%       zeros_bessel - A vector containing the first n zeros.

    zeros_bessel = zeros(n, 1);
    % Initial guess for the first zero
    x0 = 1; 
    for i = 1:n
        % fzero finds the zero of the function handle near the guess x0
        zeros_bessel(i) = fzero(@(x) besselj(order, x), x0);
        % Update the guess for the next zero, knowing they are roughly pi apart
        x0 = zeros_bessel(i) + pi; 
    end
end