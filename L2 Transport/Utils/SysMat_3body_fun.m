 function A = SysMat_3body_fun(mu,x, y, z)
%SYSMAT_3BODY_FUN(mu,x, y, z) Computes the system matrix for a rotating 3
%body system about a nominal trajectory
% Purpose: Computes the state transition matrix for the differential
% correction procedures as described by Barden SM (page 16). Matrix is of
% the form A(t) = [0, I; U'_xx 2 Omega]
%
% Inputs:
% mu        [1x1] Gravitational parameter (mass ratio parameter) M2/(M1+M2)
%           given that M1 > M2
% x         Non dimensional x position vector (can be matrix)
% y         Non dimensional x position vector (can be matrix)
% z         Non dimensional x position vector (can be matrix)
% 
% Outputs:
% A         [6x6] System Matrix 
%
% Co-dependencies
%   None
%
% Source:
%  B. Barden, "Using Stable Manifolds to Generate Transfers in the Circular
%  Restricted Problem of Three Bodies," M.S., December 1994.
% 
% Author: Alejandro Cabrales H
% 

% formulate vectors that get repeated frequently
d = (x + mu).^2 + y.^2 + z.^2;
r = (x - (1 - mu)).^2 + y.^2 + z.^2;

%% U'xx 2nd derivative of potential function

U_xx = 1 - (1 - mu)./d.^(3/2) - mu./r.^(3/2) + ...
    3*(1 - mu)*(x + mu).^2./d.^(5/2) + 3*mu*(x - (1 - mu)).^2./r.^(5/2);

U_xy = 3*(1 - mu)*(x + mu).*y./d.^(5/2) + 3*mu*(x - (1 - mu)).*y./r.^(5/2);

U_xz = 3*(1 - mu)*(x + mu).*z./d.^(5/2) + 3*mu*(x - (1 - mu)).*z./r.^(5/2);

U_yy = 1 - (1 - mu)./d.^(3/2) - mu./r.^(3/2) + ...
    3*(1 - mu).*y.^2./d.^(5/2) + 3*mu*y.^2./r.^(5/2);

U_yz = 3*(1 - mu)*y.*z./d.^(5/2) + 3*mu*y.*z./r.^(5/2);

U_zz = -(1 - mu)./d.^(3/2) - mu./r.^(3/2) + ...
    3*(1 - mu)*z.^2./d.^(5/2) + 3*mu*z.^2./r.^(5/2);

U_star = [U_xx, U_xy, U_xz;
          U_xy, U_yy, U_yz;
          U_xz, U_yz, U_zz];

% Corioli's Term due to rotation matrix
Omega = [0 1 0; -1 0 0; 0 0 0];

A = [zeros(3), eye(3); U_star, 2*Omega];    
end

