function [U] = potential_3body_fun(mu,x,y,z)
%POTENTIAL_3BODY_FUN(mu,x,y,z) Computes the pseudo potential of a function
%
% Inputs:
% mu        [1x1] Gravitational parameter (mass ratio parameter) M2/(M1+M2)
%           given that M1 > M2
% x         Non dimensional x position vector (can be matrix)
% y         Non dimensional x position vector (can be matrix)
% z         Non dimensional x position vector (can be matrix)

% Relative position vector from primary to desired vector
d = sqrt((x + mu).^2 + y.^2 + z.^2);

% Relative position vector from secondary to desired vector
r = sqrt( (x - 1 + mu).^2  + y.^2 + z.^2 );

U = 1/2 * (x.^2 + y.^2) + (1-mu)./d + mu./r;

end

