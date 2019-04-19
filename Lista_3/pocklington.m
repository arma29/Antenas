function [P]=pocklington(r, radius, lambda)
% P = pocklington(r, radius)
% r = distance from this point to the skin of the wire of the segment
% being evalutated. r may be a vector of such distances.
% radius = radius of the wire being analyzed
% lambda = wavelength
%
% This function returns the integrand for each segment of a
% wire being evalutated using Pocklingtons Ingetral Equation. If the
% input parameter r contains a vector of distances, then P is returned as
% a vector of integrand values. If r is a scalar value then, P is returned
% as a scalar.


mu0 = 4*pi*10^-7; % permeability (in newtons/amp^2)
epsilon0 = 8.854187817*10^-12; % permittivity (in Farads/meter )
eta = sqrt(mu0/epsilon0); % intrinsic impedance (in ohms)
beta = (2*pi)/lambda; % wave number (in radians/meter)
P = exp(-j*beta*r).*((1+j*beta*r).*(2*r.^2-3*radius^2)+(beta*radius*r).^2);
P = P./(4*pi*r.^5);
