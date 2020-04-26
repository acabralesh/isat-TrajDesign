function [sc_pos_bei,sc_vel_bei] = ber_to_bei(bei_R_ber, omega_ber, sc_pos_ber, sc_vel_ber)
%[sc_pos_bei,sc_vel_bei] = BER_TO_BEI(bei_R_ber,omega_ber, sc_pos_ber, sc_vel_ber)
%Converts BER to BEI position and velocity
%   BER is the Baricenter Ecliptic Rotating frame defined with origin as
%   the baricenter of the Sun-Earth system with a rotation rate, X points
%   towards earth, z is in the direction of the angular momentum. BEI is
%   the Baricenter Ecliptic Inertial frame.
%
% Inputs
% - bei_R_ber   [3x3] Rotation matrix that describe the rotation from the
%               BER frame to BEI
% - omega_ber   [3x1] Angular velocity in the BER frame 
% - sc_pos_ber  [3x1] Position of the spacecraft in BER frame
% - sc_vel_ber  [3x1] Velocity of the spacecraft in BER frame
%
% Outputs
% - sc_pos_bei  [3x1] Position of the spacecraft in BEI frame
% - sc_vel_bei  [3x1] Velocity of the spacecraft in BEI frame
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

% change 
sc_pos_bei = bei_R_ber*sc_pos_ber;

cross_omegaber = [0 -omega_ber(3) omega_ber(2);
    omega_ber(3) 0 -omega_ber(1);
    -omega_ber(2) omega_ber(1) 0];

sc_vel_bei = bei_R_ber*(sc_vel_ber + cross_omegaber*sc_pos_ber);
end

