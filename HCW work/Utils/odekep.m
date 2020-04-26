function [xdot] = odekep(t,x,force, mu, m)
% [xdot] = ODEKEP(t,x,force, SIM, SC) wrapper function for odeXY
% 
% Purpose: Wraps the kepler acceleration to produce the state vector

% Inputs
%   t           [1x1] time vector
%   x           [2*ndimx1] state vector fo the system
%   force       [ndimx1] constant force excerted for the period of time
%   SIM         Structure for SIM wide parameters
%   - UNI.mu    [1x1] Gravitational constant of the larger body
%   SC          Structure for the spacecraft orbiting
%   - PHY.m     [1x1] Mass of the spacecraft
%
% Assumptions:
%   - Assumes force is constant for the period of integration
%   
% Co-dependencies
%   kepler_accel
%
% Source:
%  Alfriend KT. Spacecraft Formation Flying?: Dynamics, Control and
%  Navigation. Elsevier Science and Technology Books, Inc., 2010
% 
% Author: Alejandro Cabrales H
% 
xdotkep = kepler_accel(mu, m, x, force);

xdot = [x(4:6); xdotkep];

end