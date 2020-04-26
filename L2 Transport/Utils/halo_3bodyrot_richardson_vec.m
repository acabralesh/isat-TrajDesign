function [x, y, z, xd, yd, zd] = halo_3bodyrot_richardson_vec(...
    t, mu, x_l, libind, Az, phi, n)
%HALO_3BODYROT_RICHARDSON_VEC(t, mu, x_l, libind, Az, phi, n) Richardson
%3rd order vector
%
% Purpose: Obtains the time evolution state vector for richardson 3rd order
% solution
%
% Inputs:
% t         [nx1] non dimensionalized time vector to compute the richardson
%           solution about
% mu        [1x1] Gravitational parameter (mass ratio parameter) M2/(M1+M2)
%           given that M1 > M2
% x_l       [1x1] Distance from secondary about libration point. Note the
%           units must coincide with that of mu
% Az        [1x1] Non dimensional amplitude of desired halo orbit in the
%           out of plane direction (z)
% libind    [1x1] Indice for selecting either libration point 1 or 2
% phi       [1x1] Angle for solution (where the solution starts) normally
%           chose zero
% n         [1x1] For choosing which type of orbits (class 1 or 3) must be
%           set to 1 or 3
%
% Outputs: 
% x         [nx1] x position (rotating frame w.r.t libration point)
% y         [nx1] y position (rotating frame w.r.t libration point)
% z         [nx1] z position (rotating frame w.r.t libration point)
% xd        [nx1] x velocity (rotating frame w.r.t libration point)
% yd        [nx1] y velocity (rotating frame w.r.t libration point)
% zd        [nx1] z velocity (rotating frame w.r.t libration point)
%
% Co-dependencies
%  c_libration(mu, gamma_l, libind, ii) needed to compute c_# values 
%
% Source:
%  [1] Richardson, D.L. Celestial Mechanics (1980) 22: 241. "Analytic
%  Construction of Periodic Orbits about the Col;inear Points."
%  https://doi.org/10.1007/BF01229511
% 
% Author: Alejandro Cabrales H

gamma_l = x_l;

% ----- Begin Function -----

% c coefficients
c2 = c_libration(mu, gamma_l, libind, 2);
c3 = c_libration(mu, gamma_l, libind, 3);
c4  = c_libration(mu, gamma_l, libind, 4);

% eigen values
% obtain the solution to lamb^4 + (2-c_2)*lamb^2 + (1+2*c_2)*(1-c_2) = 0
eig_l = roots([1 0 (2-c2) 0 (1+2*c2)*(1-c2)]);

lamb = max(imag(eig_l));

% correction coefficient for y
k = 1/(2*lamb)*(lamb^2 + 1 + 2*c2);

% Correction for conmesurate condition on z
delta = lamb^2 - c2;

% Coefficients for Richardson 3rd order

% d1 d2 coefficients
d1 = 3*lamb^2/k*(k*(6*lamb^2 - 1) - 2*lamb);

d2 = 8*lamb^2/k*(k*(11*lamb^2 - 1) - 2*lamb);

% a21 a22 a23 a24 coefficients
a21 = 3*c3*(k^2 - 2)/(4*(1 + 2*c2));

a22 = 3*c3/(4*(1 + 2*c2));

a23 = -(3*c3*lamb/(4*k*d1))*(3*k^3*lamb - 6*k*(k - lamb) + 4);

a24 = -(3*c3*lamb/(4*k*d1))*(2 + 3*k*lamb);

% d21 d31 d32 coefficients 
d21 = -c3/(2*lamb^2);

d31 = 3/(64*lamb^2)*(4*c3*a24 + c4);

d32 = 3/(64*lamb^2)*(4*c3*(a23 - d21) + c4*(4 + k^2));

% b coefficients pt 1
b21 = -(3*c3*lamb/(2*d1))*(3*k*lamb - 4);

b22 = 3*c3*lamb/d1;

b31 = (3/(8*d2))*( 8*lamb*( 3*c3*(k*b21 - 2*a23) - c4*(2 + 3*k^2)) + ...
    (9*lamb^2 + 1 + 2*c2)*( 4*c3*(k*a23 - b21) + k*c4*(4 + k^2)) );

b32 = 1/d2*( 9*lamb*( c3*(k*b22 + d21 - 2*a24) - c4) + ...
    3/8*(9*lamb^2 + 1 + 2*c2)*( 4*c3*(k*a24 - b22) + k*c4) );

% a1 a2 coefficients
a1 = -3/2*c3*(2*a21 + a23 + 5*d21) - 3/8*c4*(12 - k^2);

a2 = 3/2*c3*(a24 - 2*a22) + 9/8*c4;

% a31 a32 coefficients
a31 = -(9*lamb/(4*d2))*( 4*c3*(k*a23 - b21) + k*c4*(4 + k^2) ) + ...
    ((9*lamb^2 + 1 - c2)/(2*d2))*( 3*c3*(2*a23 - k*b21) + c4*(2 + 3*k^2));

a32 = -1/d2*( 9*lamb/4*( 4*c3*(k*a24 - b22) + k*c4 ) + ...
    3/2*(9*lamb^2 + 1 - c2)*( c3*(k*b22 + d21 - 2*a24) - c4 ) );

% s1 s2 coefficients
s1 = (1/(2*lamb*( lamb*(1 + k^2) - 2*k)))*( 3/2*c3*( 2*a21*(k^2 - 2) - ...
    a23*(k^2 + 2) - 2*k*b21) - 3/8*c4*(3*k^4 - 8*k^2 + 8));

s2 = (1/(2*lamb*( lamb*(1 + k^2) - 2*k)))*( 3/2*c3*( 2*a22*(k^2 - 2) + ...
    a24*(k^2 + 2) + 2*k*b22 + 5*d21 ) + 3/8*c4*(12 - k^2) );

% l1 l2 coefficients
l1 = a1 + 2*lamb^2*s1;

l2 = a2 + 2*lamb^2*s2;

% Generate the amplitude 
Ax = sqrt((-l2*Az^2 - delta)/l1);

% w1 w2
w1 = 0;

w2 = s1*Ax^2 + s2*Az^2;

w = 1+ w1 + w2; 

% phi, psi
psi = phi + n*pi/2;

% Plotting function

tau = t;

tau1 = lamb*w*tau + phi;

deltan = 2 - n;

x = a21*Ax^2 + a22*Az^2 - Ax*cos(tau1) + (a23*Ax^2 - a24*Az^2)*cos(2*tau1) + ...
    (a31*Ax^3 - a32*Ax*Az^2)*cos(3*tau1);

y = k*Ax*sin(tau1) + (b21*Ax^2 - b22*Az^2)*sin(2*tau1) + ...
    (b31*Ax^3 - b32*Ax*Az^2)*sin(3*tau1);

z = deltan*Az*cos(tau1) + deltan*d21*Ax*Az*(cos(2*tau1) - 3) + ...
    deltan*(d32*Az*Ax^2 - d31*Az^3)*cos(3*tau1);

xd = w*lamb*(Ax*sin(tau1) - 2*(a23*Ax^2 - a24*Az^2)*sin(2*tau1) - ...
    3*(a31*Ax^3 - a32*Ax*Az^2)*sin(3*tau1));

yd = w*lamb*(k*Ax*cos(tau1) + 2*(b21*Ax^2 - b22*Az^2)*cos(2*tau1) + ...
    3*(b31*Ax^3 - b32*Ax*Az^2)*cos(3*tau1));

zd = w*lamb*(-deltan*Az*sin(tau1) - 2*deltan*d21*Ax*Az*sin(2*tau1) - ... 
    3*deltan*(d32*Az*Ax^2 - d31*Az^3)*sin(3*tau1));

end

