function bei_C_ber = gen_beiCber_fun(pos,vel)
%[cer_C_ber] = GENERATE_BEICBER_FUN(pos,vel)
% Computes the rotation matrix from Baricenter rotating frame to eci
%   
%
% Inputs
% - pos         [3x1] Position from primary body to secondary body 
% - vel         [3x1] Velocity from primary body to secondary body 
%
% Outputs
% - bei_C_ber   [3x3] Rotation matrix from BER to BEI
%
% Assumptions:
%   -
%   
% Co-dependencies
%   - 
%
% Source:
%
% Author: Alejandro Cabrales Hernandez

pos = pos(:);
vel = vel(:);

i_hat = pos/norm(pos);
k_hat = cross(pos, vel)/norm(cross(pos, vel));
j_hat = cross(k_hat, i_hat);

bei_C_ber = [i_hat, j_hat, k_hat];

end

