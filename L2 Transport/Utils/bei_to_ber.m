function [sc_pos_ber,sc_vel_ber] = bei_to_ber(ber_R_bei, omega_ber, sc_pos_bei, sc_vel_bei)
%[sc_pos_ber,sc_vel_ber] = BEI_TO_BER(ber_R_bei, omega_ber, sc_pos_bei, sc_vel_bei)
%Converts BER to BEI position and velocity
%   BER is the Baricenter Ecliptic Rotating frame defined with origin as
%   the baricenter of the Sun-Earth system with a rotation rate, X points
%   towards earth, z is in the direction of the angular momentum. BEI is
%   the Baricenter Ecliptic Inertial frame.
%
% Inputs
% - ber_R_bei   [3x3] Rotation matrix that describe the rotation from the
%               BEI frame to BER
% - omega_ber   [3x1] Angular velocity in the BEI frame 
% - sc_pos_bei  [3x1] Position of the spacecraft in BEI frame
% - sc_vel_bei  [3x1] Velocity of the spacecraft in BEI frame
%
% Outputs
% - sc_pos_ber  [3x1] Position of the spacecraft in BER frame
% - sc_vel_ber  [3x1] Velocity of the spacecraft in BER frame
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
sc_pos_ber = ber_R_bei*sc_pos_bei;

cross_omegaber = [0 -omega_ber(3) omega_ber(2);
    omega_ber(3) 0 -omega_ber(1);
    -omega_ber(2) omega_ber(1) 0];

sc_vel_ber = ber_R_bei*(sc_vel_bei - cross_omegaber*sc_pos_bei);

end

