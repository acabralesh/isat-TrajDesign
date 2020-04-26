function [SC2_rel_pos_LVLH,SC2_rel_vel_LVLH, SC2_rel_a_LVLH, lvlh_Q_eci] = eci_to_lvlh(mu, ...
    SC1_pos_ECI, SC1_vel_ECI, SC2_pos_ECI, SC2_vel_ECI)
%[SC2_rel_pos_LVLH,SC2_rel_vel_LVLH, SC2_rel_a_LVLH, lvlh_Q_eci] =
%eci_to_lvlh(mu, ... SC1_pos_ECI, SC1_vel_ECI, SC2_pos_ECI, SC2_vel_ECI)
%Converts ECI position and velocity to LVLH frame
%   LVLH frame is defined as follos: x-axis is directed along outward radia
%   to the chaser, z-axis is in the normal to the orbital plane. y-axis
%   completes the right-handed coordinate frame
%
%   Inputs
%   - mu            [1x1] Kepler Constant (G*M)
%   - SC1_pos_ECI   [3x1] Position vector of Leader in ECI frame
%   - SC1_vel_ECI   [3x1] Velocity vector of Leader in ECI frame
%   - SC1_pos_ECI   [3x1] Position vector of Follower in ECI frame
%   - SC1_vel_ECI   [3x1] Velocity vector of Follower in ECI frame
%
%   Outputs
%   - SC1_pos_LVLH  [3x1] Position vector of Follower in LVLH frame
%   - SC1_vel_LVLH  [3x1] Velocity vector of Follower in LVLH frame 
%   - SC1_rel_a_LVLH[3x1] Acceleration vector of SC2 relative to target in
%                   lvlh frame
%   - lvlh_Q_eci    [3x3] Rotation Matrix that converts eci to lvlh
%                   coordinates
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

% defined the unit vectors for LVLH frame
i_axs = SC1_pos_ECI./norm(SC1_pos_ECI);
k_axs = hA./norm(hA);
j_axs = cross(k_axs, i_axs);

% Transformation matrix (ECI to LVLH):
lvlh_Q_eci = [i_axs'; j_axs'; k_axs'];

% angular velocity vector
Om = hA/norm(SC1_pos_ECI).^2;

% rate of change of Omega vector
Om_dot = -2*(SC1_vel_ECI'*SC1_pos_ECI)/norm(SC1_pos_ECI)^2.*Om;

% Acceleration (in ECI) of Leader and Follower
aA = -mu*SC1_pos_ECI/norm(SC1_pos_ECI).^3;

aB = -mu*SC2_pos_ECI/norm(SC2_pos_ECI).^3;

% Relative position (in ECI frame) 
SC2_posrel_ECI = (SC2_pos_ECI-SC1_pos_ECI);

% Converting to LVLH frame
SC2_rel_pos_LVLH = lvlh_Q_eci*SC2_posrel_ECI;

% Relative velocity of the (in ECI frame)
SC2_velrel_ECI = (SC2_vel_ECI - SC1_vel_ECI - ...
    cross(Om, SC2_posrel_ECI));

% converting to LVLH frame
SC2_rel_vel_LVLH = lvlh_Q_eci*SC2_velrel_ECI;

% Relative acceleration of follower relative to target in ECI frame
a_rel = aB-aA - cross(Om_dot, SC2_posrel_ECI) - ...
    cross(Om, cross(Om, SC2_posrel_ECI)) - 2*cross(Om, SC2_velrel_ECI);

% Change the acceleration in LVLH
SC2_rel_a_LVLH = lvlh_Q_eci*a_rel;



end

