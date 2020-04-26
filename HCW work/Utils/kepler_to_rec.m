function [pos,vel] = kepler_to_rec(mu,a, eccentricity, inclination, longnode, argperi, meananom)
%[pos,vel] = KEPLER_TO_REC(mu,a, eccentricity, inclination, longnode,
%argperi, meananom) Converts keplerian orbital elements to
%rectangular
%   Inputs:
%   - mu            [1x1] Keplerian constant (G*m)
%   - a             [1x1] semi-major axis
%   - eccentricity  [1x1] eccentricity
%   - inclination   [1x1] inclination
%   - longnode      [1x1] RAAN
%   - argperi       [1x1] argument of perigee
%   - meananom      [1x1] mean anomaly
% 
%   Output:
%   - pos           [3x1] Position vector 
%   - vel           [3x1] Velocity vector
%
% Sources:
% J. Wisdom, Planetary Dynamical Systems. To be published. Ch. 3 Kepler
% Problem

% Solve's Kepler Equation
E = KeplerEquation(eccentricity, meananom, 1e-7);
cosE = cos(E);
sinE = sin(E);

foo = sqrt(1 - eccentricity*eccentricity);
meanmotion = sqrt(mu/(a.^3));
x = a*(cosE-eccentricity);
y = foo*a*sinE;
z = 0.0;

vx = -a*meanmotion*sinE/(1.0 - eccentricity*cosE);
vy = foo*a*meanmotion*cosE/(1.0 - eccentricity*cosE);
vz = 0.0;

pos = [x;y;z];
vel = [vx;vy;vz];

% perform the rotation transformations
pos1 = rotate_z(cos(argperi), sin(argperi), pos);
vel1 = rotate_z(cos(argperi), sin(argperi), vel);

pos2 = rotate_x(cos(inclination), sin(inclination), pos1);
vel2 = rotate_x(cos(inclination), sin(inclination), vel1);

pos = rotate_z(cos(longnode), sin(longnode), pos2);
vel = rotate_z(cos(longnode), sin(longnode), vel2);



end


% function E = keplerEq(M,e,eps)
% % Function solves Kepler's equation M = E-e*sin(E)
% % Input - Mean anomaly M [rad] , Eccentricity e and Epsilon 
% % Output  eccentric anomaly E [rad]. 
%    	En  = M;
% 	Ens = En - (En-e*sin(En)- M)/(1 - e*cos(En));
% 	while ( abs(Ens-En) > eps )
% 		En = Ens;
% 		Ens = En - (En - e*sin(En) - M)/(1 - e*cos(En));
%     end
% 	E = Ens;
% end

function E = KeplerEquation(e,M, eps)
% Solves the Kepler equation

E = M;
E = M + e*sin(E);
E = M + e*sin(E);
error = 1;

while error > eps
   esinE = e*sin(E);
   ecosE = e*cos(E);
   f = E - (M + esinE);
   fp = 1 - ecosE;
   fpp = esinE;
   fppp = ecosE;
   DeltaNewton = -f/fp;
   DeltaHalley = -f/(fp + DeltaNewton*fpp/2.0);
   DeltaDanby = -f/(fp + DeltaHalley*(fpp + DeltaHalley*fppp/3.0)/2.0);
   E = E + DeltaDanby;
   error = abs(DeltaDanby);
end

end

function Rvector = rotate_z(cosangle, sinangle, vector)
Rvector = zeros(3,1);
Rvector(1) = vector(1)*cosangle - vector(2)*sinangle;
Rvector(2) = vector(1)*sinangle + vector(2)*cosangle;
Rvector(3) = vector(3);

end

function Rvector = rotate_x(cosangle, sinangle, vector)
Rvector = zeros(3,1);

Rvector(1) = vector(1);
Rvector(2) = vector(2)*cosangle - vector(3)*sinangle;
Rvector(3) = vector(2)*sinangle + vector(3)*cosangle;

end