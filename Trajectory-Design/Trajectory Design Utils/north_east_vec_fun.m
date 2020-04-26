function [dn, de] = north_east_vec_fun(pos_launch_vec)
% NORTH_EAST_vec_fun obtains the north and earth unit vectors from position
% of earth in ECEF coordiantes
%
% Inputs:
% pos_launch_vec    [3x1] vectors that describes the position of launch in
%                   ECEF coordinates
% 
% Outputs:
% dn                [3x1] Unit vector representing the north direction on a
%                   sphere
% de                [3x1] Unit vector representing the east direction on a
%                   sphere
%
% Co-dependencies
%   None
%
% Source:
%  https://www.movable-type.co.uk/scripts/latlong-vectors.html
% 
% Author: Alejandro Cabrales H
% 

% north pole in ECEF
n_ecef = [0 0 1]'; % nort h pole

a_vec = pos_launch_vec/norm(pos_launch_vec);

de = cross(n_ecef, a_vec); % tangential due east vector

% due east vector
de = de/norm(de);

% due north vector  from launch
dn = cross(a_vec, de);

dn = dn/norm(dn);

de = de(:);
dn = dn(:);
end
