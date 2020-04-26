function [x_l1,x_l2] = find_librationpts(mu)
%[x_l1,x_l2] = FIND_LIBRATIONPTS(mu) 
%
% Purpose: Returns the location of both of the libration points
%
% Inputs:
% mu        [1x1] Gravitational parameter (mass ratio parameter) M2/(M1+M2)
%           given that M1 > M2
%
% Outputs: 
% x_l1      [1x1] distance of L1 relativive to secondary orbit
% x_l2      [1x1] distance of L2 relativive to secondary orbit
%
% Co-dependencies
%  None
%
% Source:
%  [1] Richardson, D.L. Celestial Mechanics (1980) 22: 241. "Analytic
%  Construction of Periodic Orbits about the Col;inear Points."
%  https://doi.org/10.1007/BF01229511
% 
% Author: Alejandro Cabrales H

% Find Location of L1 L2 based on roots of polinomial 
L1 = roots([1,-(3-mu),(3-2*mu),-mu,2*mu,-mu]);
L2 = roots([1,(3-mu),(3-2*mu),-mu,-2*mu,-mu]);

% x-position of Lagrange points in rotating frame
% Note xL1 is non dimensional
x_l1 = L1(imag(L1)==0 & real(L1)>0);
x_l2 = L2(imag(L2)==0 & real(L2)>0);

end

