function [sc_pos_cei,sc_vel_cei] = bei_to_cei(cei_R_bei, omega_ber,r_earth_bei, sc_pos_bei, sc_vel_bei)
%[sc_pos_cei,sc_vel_cei] = bei_to_cei(cei_R_bei, omega_ber,r_earth_bei,
%sc_pos_bei, sc_vel_bei) Converts BEI to CEI position and velocity
%   BEI is the Baricenter Ecliptic Inertial frame defined with origin as
%   the baricenter X points towards Aires, Z points towards the angular
%   momentum of the Earth Sun system, and Y completes the axis. CEI is the
%   same just centered at Earth with the velocity of earth taken out
%
% Inputs
% - cei_R_bei   [3x3] Rotation matrix that describe the rotation from the
%               BEI frame to CEI
% - omega_ber   [3x1] Angular velocity in the BEI frame 
% - r_earth_bei [3x1] position of earth in BEI frame
% - sc_pos_bei  [3x1] Position of the spacecraft in BEI frame
% - sc_vel_bei  [3x1] Velocity of the spacecraft in BEI frame
% - r_earth     [3x1] position of earth
%
% Outputs
% - sc_pos_cei  [3x1] Position of the spacecraft in CEI frame
% - sc_vel_cei  [3x1] Velocity of the spacecraft in CEI frame
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

sc_pos_cei = cei_R_bei*(sc_pos_bei - r_earth_bei);

cross_omegaber = [0 -omega_ber(3) omega_ber(2);
    omega_ber(3) 0 -omega_ber(1);
    -omega_ber(2) omega_ber(1) 0];

sc_vel_cei = cei_R_bei*(sc_vel_bei - cross_omegaber*r_earth_bei);
