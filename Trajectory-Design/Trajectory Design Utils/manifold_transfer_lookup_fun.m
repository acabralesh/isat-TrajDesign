function [delta_vtot_min,trans_duration_min, mani_duration_min, out_info] = manifold_transfer_lookup_fun(launch_state,manifolds, bei_C_ber, convert)
%[delta_v,trans_duration, mani_duration] = MANIFOLD_TRANSFER_LOOKUP_FUN(launch_state,manifolds, bei_C_ber, convert)
%Obtains the cheapest trajectory to a manifold of the desired body
%   
% Inputs
%  - launch_state       [6x1] Launch state vector (1:3 = position) (4:6 =
%                       velocity. Units must be in SI units
%  - manifolds          [struct] Structure with manifolds values
%    - sol              struct] with individual manifold trajectory
%      - xnumvec        [nx6] Vector for manifold trajectory. NOTE: initial
%                       position is on HALO orbit. Untis are dimensionless
%     - t               [nx1] Time vector (time runs negative). Unitless
%  - mu_value           [1x1] mu value for system in question
%  - bei_C_ber          [3x3] Rotation matrix from Barycenter rotating
%                       frame to inertial frame based on JD vector and
%                       planet ephemeris
%  - convert            [struct] Matlab structure that converts time
%                       position and velocity to SI units
%
% Outputs
%  - delta_vtot_min     [1x1] delta-v cost for the mission in SI untis
%  - trans_duration_min [1x1] Time for Lamberts transfer in seconds
%  - mani_duration_min  [1x1] Time to reach manifold departure point
%  - out_info           [struct] structure with information about transfer
%
% Assumptions:
% -NOTE: The current assumption is that the lamber solution is correct.
% This works for a 2 body system, however, due to the 3rd body effects this
% is not the case. Future work will involve using differential corrector
% schemes
%   
% Co-dependencies
% - lambert_mina    Computes the solution to lamberts problem
% - ber_to_bei      Transform from Rotating frame to Inertial frame
% - bei_to_cei      Transform from Barycenter to Earth Center frame
%
% Source:
%
% Author: Alejandro Cabrales Hernandez


% Need to convert positions

%% Unpack
mu = manifolds.mu;
r_earth_ber = manifolds.earth_ber;

pos_units_convert   = convert.pos;
vel_units_convert   = convert.vel;
time_units_convert  = convert.time;

Gmu = 3.986004418*10^(14)*METERS^3/SECONDS^2; % Sorry for the hardcode

%% Values needed for converting 
omega_n = [0 0 1]';



%% Set up the for loop values
delta_vtot_min = Inf;
delta_v1_min = Inf;
delta_v2_min = Inf;
v1_min = Inf;
v2_min = Inf;
trans_duration_min = NaN;
mani_duration_min = NaN;
target_state_min = NaN*zeros(6,1);
mani_id_min = NaN;


% Run over each manifolds
for ii = 1:(14+0*length(manifolds.sol))
    
    xnumvec     = manifolds.sol{ii}.xnumvec;
    tnumvec     = -manifolds.sol{ii}.tnumvec;
    
    target_state = reshape(xnumvec(end,:),6,1);
    target_time = tnumvec(end);
    
    % TODO - Make sure frame transformations are done correctly.
    % Convert to Inertial frame
    [target_pos_bei,target_vel_bei] = ber_to_bei(bei_C_ber, omega_n,target_state(1:3), target_state(4:6));
    
    % Obtain Location of Earth:    
    r_earth_bei = bei_C_ber*r_earth_ber;
    
    % Convert to Earth Centered frame (Translation + velocity of
    % earth)
    [target_pos_eci,target_vel_eci] = bei_to_cei(eye(3), bei_C_ber*omega_n,r_earth_bei, target_pos_bei, target_vel_bei);

    % Convert to correct Units:
    target_state_eci = [target_pos_eci*pos_units_convert; target_vel_eci*vel_units_convert];
    
    
    % ---- Lamberts -----
    [v1_lambvec,v2_lambvec, T_lamb] = lambert_mina(Gmu, launch_state, target_state_eci);
    
    % TODO in future: Use Differential Corrector scheme to guarantee that
    % you reach desired state
    
    % Check for colision:
    % Collision to earth will only happen if the projection of the
    % deltav1_lamvec onto launch_state.position is negative
    
    col = sign(dot(v1_lambvec, launch_state(1:3)));
    
    % the results from lambert_mina are the velocity of the transfer
    % trajectory, we are interested in computing the deltav1 velocitites
    deltav1_lambvec = v1_lambvec-launch_state(4:6);
    deltav2_lambvec = target_state_eci(4:6) - v2_lambvec;
    
    
    delta_v_inloop = norm(deltav1_lambvec) + norm(deltav2_lambvec);
    
    %  if there is no colision and the delta-v is much larger then use thsi
    %  method
    if col > 0 && delta_v_inloop <= delta_vtot_min
        
        % update the current best estimates
        delta_vtot_min      = delta_v_inloop;
        v1_min              = v1_lambvec;
        v2_min              = v2_lambvec;
        delta_v1_min        = deltav1_lambvec;
        delta_v2_min        = deltav2_lambvec;
        trans_duration_min  = T_lamb;
        mani_duration_min   = target_time*time_units_convert;
        mani_id_min         = ii;
        target_state_min    = target_state_eci;
        
    end
    
    
end

% save all data forward
out_info.delta_vtot_min = delta_vtot_min;
out_info.delta_v1_min = delta_v1_min;
out_info.delta_v2_min = delta_v2_min;
out_info.v1_min = v1_min;
out_info.v2_min = v2_min;
out_info.trans_duration_min = trans_duration_min;
out_info.mani_duration_min = mani_duration_min;
out_info.mani_id_min = mani_id_min;
out_info.target_state_min = target_state_min;

end
