function [value,isterminal,direction] = diffcorr_3body_zx_cross_fun(t,xvec)
%DIFFCORR_ZX_CROSS_FUN(t,xvec)Event location for differential correction
%
% Purpose: Sets the event function for the ODEXY in which it commands ode
% to stop evaluation once the orbit crosses the xz plane
%
% Inputs:
% t         [1x1] nondimensionalized time
% xvec      [42x1] nondimensionalized state vector + flattened STM vector.
%           First 6 elements corresponds to the state, next 36 elements are
%           the elements corresponding to the state transition matrix. For
%           example if STM = magic(3), the flattened vectors are STM(:)
%           which corresponds to flattening the vector column wise
%
% Outputs: 
% value         [1x1] Selects the value for which the event
%               corresponds (y = 0) 
% isterminal    [1x1] sets equal to 1, indicates ode45 will terminate
%               evaluation
% direction     [1x1] for now, direction doesnt matter (useful for Class 1
%               Class 2 orbits)
%
% Co-dependencies
%  None
%
% Source:
%  https://www.mathworks.com/help/matlab/math/ode-event-location.html
% 
% Author: Alejandro Cabrales H
% 
value = xvec(2) ; %t - 1.50623074975922;     % Detect y = 0
isterminal = 1;      % Stop the integration
direction = 0;       % No direction matters
end

