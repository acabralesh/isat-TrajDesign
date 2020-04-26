function xdotvec = eom_3bodystate_STM_fun(t,xvec, mu)
%EOM_3BODYSTATE_STM_FUN(t,xvec, mu) Outputs the vector xdot = f(x)for use
%in integration
% Purpose: Computes the non linear equations of motion of the combined
% state vector and state transition matrix for the 3rd body problem in
% baricenter rotating coordinates
%
% Inputs:
% t         [1x1] nondimensionalized time
% xvec      [42x1] nondimensionalized state vector + flattened STM vector.
%           First 6 elements corresponds to the state, next 36 elements are
%           the elements corresponding to the state transition matrix. For
%           example if STM = magic(3), the flattened vectors are STM(:)
%           which corresponds to flattening the vector column wise
% mu        [1x1] Gravitational parameter (mass ratio parameter) M2/(M1+M2)
%           given that M1 > M2
%
% Outputs:
% xdotvec   [42x1] State + STM vector. Derivative of xvec
%
% Co-dependencies
%   eom_3_body_rot_fun  Function that returns the first 6 elements of the
%                       the xdotvec (computes the equation of motion for
%                       rotating 3rd body problem)
%   SysMat_3body_fun    Computes the matrix Phidot = A(t) phi which helps
%                       in obtaining the other 36 elements
%
% Source:
%  [1] B. Barden, "Using Stable Manifolds to Generate Transfers in the
%  Circular Restricted Problem of Three Bodies," M.S., December 1994.
%  [2] Richardson, D.L. Celestial Mechanics (1980) 22: 241. "Analytic
%  Construction of Periodic Orbits about the Col;inear Points."
%  https://doi.org/10.1007/BF01229511
% 
% Author: Alejandro Cabrales H
% 

% generate stateeom

xstateeom = eom_3body_rot_fun(t, xvec(1:6), mu);

% compute the system matrix phi
A = SysMat_3body_fun(mu, xvec(1), xvec(2), xvec(3));

% regenerate the phi vector
phi = reshape(xvec(7:end), 6,6); 

% compute  phidot
STMdot = A*phi;

xSTMeom = STMdot(:);

xdotvec = [xstateeom ; xSTMeom];

end

