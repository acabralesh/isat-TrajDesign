function c = c_libration(mu, gamma_l, libind, n)
%C_LIBRATION(mu, gamma_l, libind, n) Computes the c(n) value for use in
%Libration approximation
% Purpose:  Computes the coefficient c resulting from linear approximation
% of the 3rd body problem in search for libration points
% 
% Inputs:
% mu        [1x1] Gravitational parameter (mass ratio parameter) M2/(M1+M2)
%           given that M1 > M2
% gamma_l   [1x1] dimensionless quantity expressing distance from secondary
%           (earth) to the location of lagrange point
% libind    [1x] Index corresponding to wether we are interested in
%           libration point 1 or 2
% n         [1x1] nth term of the c expansion
%
% Outputs:
% c         [1x1] nth term in c expansion
%
% Co-dependencies
%   None
%
% Source: Richardson, D.L. Celestial Mechanics (1980) 22: 241. "Analytic
% Construction of Periodic Orbits about the Col;inear Points."
% https://doi.org/10.1007/BF01229511
% 
% Author: Alejandro Cabrales H
%  

% decision variable whehter using L1 or L2
lib = [1, -1];

libcorr = lib(libind);

c = 1./gamma_l.^3 *( (libcorr).^n*mu + ...
    (-1).^n*(1-mu).*gamma_l.^(n+1)./(1 -1*libcorr*gamma_l).^(n+1));


end

