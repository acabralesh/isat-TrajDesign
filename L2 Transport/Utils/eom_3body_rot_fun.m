function [xdotvec] = eom_3body_rot_fun(t, xvec, mu)
%EOM_3BODY_ROT_FUN(t, xvec, mu) equations of motion for 3 body problem,
%normalized
%
% Inputs:
% t         [1x1] time vector
% xvec      [6x1] Non dimensional state (can be matrix)
% mu        [1x1] Gravitational parameter (mass ratio parameter) M2/(M1+M2)
%           given that M1 > M2

% Unpack values
x = xvec(1);
y = xvec(2);
z = xvec(3);
xd = xvec(4);
yd = xvec(5);
zd = xvec(6);

% either use diff_potential (slower) or manual
% dudxyz = diff_potential_3body_fun(mu, [x;y;z], 1:3);

% xdotvec = [xd; yd; zd; dudxyz(1) + 2*yd; dudxyz(2) - 2*xd; dudxyz(3)];


% manual
% Relative position vector from primary to desired vector
d = sqrt((x + mu).^2 + y.^2 + z.^2);

% Relative position vector from secondary to desired vector
r = sqrt((x - (1 - mu)).^2  + y.^2 + z.^2 );

dudx = x - (1 - mu).*(x + mu)./d.^3 - mu*(x - (1 - mu))./r.^3;
dudy = y - (1 - mu).*y./d.^3 - mu*y./r.^3;
dudz = -(1 - mu)*z./d.^3 - mu*z./r.^3;

xdotvec = [xd; yd; zd; dudx + 2*yd; dudy - 2*xd; dudz];


            
end

