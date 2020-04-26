function Phi = stm_HCW(n, t)
%Phi = stm_HCW(n, t) Computes the state transition matrix for HCW dynamics
%
% Purpose: Generates the 6x6 matrix state transition matrix Phi(t,t0) where
% t0 is assumed to be 0
%
% Inputs:
% n         [1x1] Mean motion of leader n = sqrt(mu/a^3)
% t         [1x1] time 
% 
% Outputs:
% Phi       [6x6] State transition matrix
%
% Co-dependencies
%   None
%
% Source:
%  Alfried K., Vadali S.. "Spacecraft Formation Flying: Dynamics Control
%  and Navigation." First Edition. Butterworth-Heinemann. 2009. Pg. 84-86.
% 
% Author: Alejandro Cabrales H
% 

cn = cos(n*t); sn = sin(n*t);

Phi = [4-3*cn, 0, 0, sn/n, 2/n-2*cn/n, 0;
    -6*n*t + 6*sn, 1, 0, -2/n + 2*cn/n, 4*sn/n - 3*t, 0;
    0, 0, cn, 0, 0, sn/n;
    3*n*sn, 0, 0, cn, 2*sn, 0;
    -6*n + 6*n*cn, 0, 0, -2*sn, -3 + 4*cn, 0;
    0, 0, -n*sn, 0, 0, cn];


    


end

