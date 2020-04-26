function [xdotdot] = kepler_accel(mu,m, x, force)
% kepler_force: Computes the acceleration due to the kepler 2 body problem
% Inputs
%   SIM         Structure for SIM wide parameters
%   - mu        Gravitational constant of the larger body
%   SC          Structure for the spacecraft orbiting
%   - m         Mass of the spacecraft
%   x           (3x1 vector) State pos. vector 
%   force       (3x1 vecotr) constant force excerted for the period of time
% 
% Assumptions:
%   - Assumes that the gravitating bodies are spherical
%   - There are no tidal forces
%   - The primary mass is much larger than orbiting body mass
%   
% Co-dependencies
%   None
%
% Source:
%  Alfriend KT. Spacecraft Formation Flying?: Dynamics, Control and
%  Navigation. Elsevier Science and Technology Books, Inc., 2010
% 
% Author: Alejandro Cabrales H
% 

% --- Unpack values ---
% mu = SIM.UNI.mu;
% m = SC.PHY.mass;

% Compute the position vector magnitude
r = sqrt(x(1).^2 + x(2).^2 + x(3).^2);
 
% xdotdot(1) = -mu*x(1)./r.^3 + force(1)/m;
% xdotdot(2) = -mu*x(2)./r.^3 + force(2)/m;
% xdotdot(3) = -mu*x(3)./r.^3 + force(3)/m;

xdotdot = -mu*x(1:3)./r.^3 + force/m;

% xdotdot = xdotdot';

end

