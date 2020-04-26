function [sc_pos_bei,sc_vel_bei] = cei_to_bei(bei_R_cei, omega_ber,r_earth_bei, sc_pos_cei, sc_vel_cei)
%[sc_pos_bei,sc_vel_bei] = CEI_TO_BEI(bei_R_cei, omega_ber,r_earth_bei, sc_pos_cei, sc_vel_cei)
%Converts CEI to BEI position and velocity
%   BEI is the Baricenter Ecliptic Inertial frame defined with origin as
%   the baricenter X points towards Aires, Z points towards the angular
%   momentum of the Earth Sun system, and Y completes the axis. CEI is the
%   same just centered at Earth with the velocity of earth taken out
%
% Inputs
% - bei_R_cei   [3x3] Rotation matrix that describe the rotation from the
%               CEI frame to BEI
% - omega_ber   [3x1] Angular velocity in the BEI frame 
% - r_earth_bei [3x1] Position of Earth in BEI frame
% - sc_pos_cei  [3x1] Position of the spacecraft in CEI frame
% - sc_vel_cei  [3x1] Velocity of the spacecraft in CEI frame
%
% Outputs
% - sc_pos_bei  [3x1] Position of the spacecraft in CEI frame
% - sc_vel_bei  [3x1] Velocity of the spacecraft in CEI frame
%
% Assumptions:
%   - 
%   
% Co-dependencies
%   None
%
% Source:
%
% Author: Alejandro Cabrales Hernandez

sc_pos_bei = bei_R_cei*sc_pos_cei + r_earth_bei;

cross_omegaber = [0 -omega_ber(3) omega_ber(2);
    omega_ber(3) 0 -omega_ber(1);
    -omega_ber(2) omega_ber(1) 0];

sc_vel_bei = bei_R_cei*sc_vel_cei + cross_omegaber*r_earth_bei;

end