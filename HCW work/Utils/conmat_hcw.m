function [A,B] = conmat_hcw(n, dim)
%[A,B] = conmat_hcw(n) Computes the continuous system matrices A, B
%
% Purpose: Generates the system matrix for the hcw problem in 2 or 3d
%
% Inputs:
% n         [1x1] Mean motion of leader n = sqrt(mu/a^3)
% dim       [1x1] either 2 or 3 representing 2D or 3D motion
% 
% Outputs:
% A       [2dimx2dim] A matrix for xdot = A x + B u system
% B       [2dimxdim] B matrix for same system
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

switch dim
    
    case 2
        
        B = [0 0; 0 0; 1 0; 0 1];
        
        A = [zeros(2,2), eye(2,2); 3*n^2 0 0 2*n; 0 0 -2*n 0];
        
    case 3
        
        B = [zeros(3,3); eye(3,3)];
        
        A = [zeros(3,3), eye(3,3);
            3*n^2 0 0 0 2*n 0;
            0 0 0 -2*n 0 0;
            0 0 -n^2 0 0 0];
end

end

