function [a, eccentricity, inclination, longnode, argperi, meananom] = rec_to_kepler(mu, pos, vel)
%	[a, eccentricity, inclination, longnode, argperi, meananom] =
%	rec_to_kepler(mu, pos, vel) Converts rectangular coordinates to
%	keplerian elements
%   Inputs:
%   - mu            [1x1] Keplerian constant (G*m)
%   - pos           [3x1] Position vector 
%   - vel           [3x1] Velocity vector
% 
%   Output:
%   - a             [1x1] semi-major axis
%   - eccentricity  [1x1] eccentricity
%   - inclination   [1x1] inclination
%   - longnode      [1x1] RAAN
%   - argperi       [1x1] argument of perigee
%   - meananom      [1x1] mean anomaly
%
%
% Sources:
% J. Wisdom, Planetary Dynamical Systems. To be published. Ch. 3 Kepler
% Problem

% Computes the cross product 
r_cross_v = cross(pos,vel);

% sum of squares of angular momentum
hs = sum(r_cross_v.*r_cross_v); 

% magnitude of angular momentum
h = sqrt(hs);

% magnitude of position vector
r = norm(pos);

% sum of squares for velocity
vs = sum(vel.*vel);

% dot product of position to velocity (used for . . . )
rdotv = dot(pos, vel);

rdot = rdotv/r;

param = hs/mu;

cosi = r_cross_v(3)/h;

inclination = acos(cosi);

% Computes the longitude node
longnode = wrapTo2Pi(atan2(r_cross_v(1), -r_cross_v(2)));

a = 1/(2/r - vs/mu);

% e*trueanomaly 
ecostrueanom = param/r - 1;

esintrueanom = rdot*h/mu;

% obtain the eccentricity of the orbit
eccentricity = sqrt(ecostrueanom.^2 + esintrueanom.^2); 

trueanom = wrapToPi(atan2(esintrueanom, ecostrueanom));

cosnode = cos(longnode);

sinnode = sin(longnode);

% checks for when cosine of inclination is too small
if abs(cosi) > 0.001
   rcosu = pos(1)*cosnode + pos(2)*sinnode;
   rsinu = (pos(2)*cosnode - pos(1)*sinnode)/cosi;
else
    deltai = 1;
    state2 = rotate_z(cos(longnode), -sin(longnode), pos);
    state3 = rotate_x(cos(deltai), sin(deltai), state2);
    pos = rotate_z(cos(longnode), sin(longnode), state3);
    rcosu = pos(1)*cosnode + pos(2)*sinnode;
    rsinu = (pos(2)*cosnode - pos(1)*sinnode)/cos(inclination+deltai);
end

argperi = wrapTo2Pi(atan2(rsinu, rcosu)-trueanom);
eccanom = 2*atan(sqrt((1 - eccentricity)/(1 + eccentricity))*tan(trueanom/2.0));
meananom = wrapTo2Pi(eccanom - eccentricity*sin(eccanom));
end

function Rvector = rotate_z(cosangle, sinangle, vector)
Rvector = zeros(3,1);
Rvector(1) = vector(1)*cosangle - vector(2)*sinangle;
Rvector(2) = vector(1)*sinagle + vector(2)*cosangle;
Rvector(3) = vector(3);

end

function Rvector = rotate_x(cosangle, sinangle, vector)
Rvector = zeros(3,1);

Rvector(1) = vector(1);
Rvector(2) = vector(2)*cosangle - vector(3)*sinangle;
Rvector(3) = vector(2)*sinagle + vector(3)*cosangle;

end