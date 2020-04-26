function [SC2_pos_ECI,SC2_vel_ECI, eci_Q_lvlh] = lvlh_to_eci(SC1_pos_ECI, ...
    SC1_vel_ECI,SC2_rel_pos_LVLH,SC2_rel_vel_LVLH)
%[SC2_pos_ECI,SC2_vel_ECI, eci_Q_lvlh] = LVLH_TO_ECI(SC1_pos_ECI,
%SC1_vel_ECI,SC2_rel_pos_LVLH,SC2_rel_vel_LVLH) Converts relative position
%in LVLH frame to ECI frame
%   LVLH frame is defined as follos: x-axis is directed along outward radia
%   to the chaser, z-axis is in the normal to the orbital plane. y-axis
%   completes the right-handed coordinate frame
%
%   Inputs
%   - mu            [1x1] Kepler Constant (G*M)
%   - SC1_pos_ECI   [3x1] Position vector of Leader in ECI frame
%   - SC1_vel_ECI   [3x1] Velocity vector of Leader in ECI frame
%   - SC1_pos_LVLH  [3x1] Position vector of Follower in LVLH frame
%   - SC1_vel_LVLH  [3x1] Velocity vector of Follower in LVLH frame 
%
%   Outputs
%   - SC2_pos_ECI   [3x1] Position vector of Follower in ECI frame
%   - SC2_vel_ECI   [3x1] Velocity vector of Follower in ECI frame
%
% Assumptions:
%   - Orbit should be circular
%   
% Co-dependencies
%   None
%
% Source:
%   CURTIS, H. D. Orbital mechanics for engineering students. [S.l.]:
%   Oxford?: Butterworth-Heinemann, 2013., 2013. Ch 7.

% Angular momentum vector for Leader orbit
hA = cross(SC1_pos_ECI, SC1_vel_ECI);

% angular velocity vector
Om = hA/norm(SC1_pos_ECI).^2;

% defined the unit vectors for LVLH frame
i_axs = SC1_pos_ECI./norm(SC1_pos_ECI);
k_axs = hA./norm(hA);
j_axs = cross(k_axs, i_axs);

% Transformation matrix (ECI to LVLH):
lvlh_Q_eci = [i_axs'; j_axs'; k_axs'];

% transformation matrix (LVLH to ECI);
eci_Q_lvlh = lvlh_Q_eci';

SC2_pos_ECI = eci_Q_lvlh*SC2_rel_pos_LVLH + SC1_pos_ECI;

SC2_vel_ECI = eci_Q_lvlh*SC2_rel_vel_LVLH + SC1_vel_ECI + ...
    cross(Om, eci_Q_lvlh*SC2_rel_pos_LVLH);


end

