function [delta_v, t_launch, t_duration, out_struct] = traj_lookup(t_start, t_end, target_loc, input_struct)
%[delta_v, t_launch, duration] = TRAJ_LOOKUP(t_start, t_end, destination, input_struct)
%Obtains the trajectory Delta-V, launch time, and travel duration to get to
%the desired destination
%   
%
% Inputs
% - t_start     [1x3] Time vector in [Year, Month, Day]for beginning of 
%               lookup time for the optimization run 
% - t_end       [1x3] Time vector in [Year, Month, Day] for end of lookup
%               time for the optimization run 
% - target_loc  [str] String for either 'EML1' or 'SEL2' as the desired
%               destination
% - input_struct  [TBD] Input flag for fancy stuff
%
% Outputs
% - delta_v     [1x1] Total cost of delta-v trip to achieve the target
% - t_launch    [1x6] Time vector [Year, Moth, Day, Hour, Minute, Second]
% - t_duration  [1x1] Time duration to reach destination in days
% - out_struct  [struct] output struct with fancy stuff
%
% Assumptions:
%   - Launch solution only looks at discrete number of trajectories.
%       Impact: Solutions are not true optimal and there is almost
%       certainty that a better trajectory is available
%   - No Solar Radiation Pressure or 4th body effects
%       Impact: Trajectory currently does not have Solar Radiation pressure
%       which will require a mid-course correction to maintain nominal
%       orbit. Furthermore, the 4th body effect (moon or sun) is not
%       present
%   - No differential corrector in place (This is the next focus)
%       Impact: true delta-v to reach desire manifold might change (change
%       is <5% based on some testing scenarios
%   
% Co-dependencies
%   - UNITS folder is in path
%   - juliandate -- (Aerospace Blockset) for conputing time to julian date
%   - planetEphemeris -- (aeroDataPackage) Download required
%   - ECEFtoECI -- (Matlab package) Download Require
%   - gen_beiCber_fun -- Located in Utils 
%   - manifold_transfer_lookup -- Located in Utils has several dependencies
%
% Source:
%
% Author: Alejandro Cabrales Hernandez

%% Unpack input_flag

% TODO use varagin
% Obtain Launch position
lat = input_struct.launch_lat;
long = input_struct.launch_long;
launch_pos_ecef =lla2ecef( [lat/DEGREES, long/DEGREES, 0]);

% Unpack Launch Atmospheric Drag
atm_drag = input_struct.atmosphere_drag; 

% Print Values
verbose = input_struct.verbose;
%% Set Up Units conversions

% Compute the Julian traget % TODO double check the input vector

JD_start    = juliandate(t_start);
JD_end      = juliandate(t_end); 

switch target_loc
    
    case 'SEL2'
        
        if verbose
            fprintf('Obtaining Manifolds to SEL2 ...')
        end
        
        % load the SEL2 Manifolds
        load('SEL2_manifolds.mat');
        
        % Set primary, secondary bodies
        primary_body    = 'Sun';
        secondary_body  = 'Earth';
        
        % --- Unpack values for SEL2 ---
        manifolds.a         = SEL2_manifolds.a;
        manifolds.mu        = SEL2_manifolds.mu;
        manifolds.sol       = SEL2_manifolds.sol;
        manifolds.earth_ber = [1-SEL2_manifolds.mu; 0; 0];
        
        %  --- Unit conversions ----
        % This converts the nondimensional value to dimensional value        
        pos_units_convert   = manifolds.a;                  %[m]
        vel_units_convert   = manifolds.a/(1*YEARS)*2*pi;    %[m/s]
        time_units_convert  = (1*YEARS)/(2*pi);          %[s]
        
        if verbose
            fprintf('done\n')
        end
    case 'EML1'
        
        if verbose
            fprintf('Obtaining Manifolds to SEL2 ...')
        end
        
        % load the EML1 Manifolds
        load('EML1_manifolds.mat');
        
        
        % Set primary, secondary bodies
        primary_body = 'Earth';
        secondary_body = 'Moon';

        % --- Unpack values for SEL2 ---
        manifolds.a         = EML1_manifolds.a;
        manifolds.mu        = EML1_manifolds.mu;
        manifolds.sol       = EML1_manifolds.sol;
        manifolds.earth_ber = [-EML1_manifolds.mu; 0; 0];
        
        %  --- Unit conversions ----
        % This converts the nondimensional value to dimensional value
        
        pos_units_convert   = manifolds.a ;                       %[m]
        vel_units_convert   = manifolds.a /(27.345*DAYS)*2*pi;    %[m/s]
        time_units_convert  = (27.345*DAYS)/(2*pi);               %[s]
        
        if verbose
            fprintf('done\n')
        end
    otherwise
        error('Currently only SEL2 and EML1 system are supported');
end

% Store conversion data
convert.pos         = pos_units_convert;
convert.vel         = vel_units_convert;
convert.time        = time_units_convert;

%% Set up Optimization loop
delta_v_best        = Inf;
t_launch_best       = NaN;
t_duration_best     = NaN;
out_info_best       = NaN;
% Compute the time vector for which the optimization loop will run
time_opt_vec        = JD_start:1*HOURS/DAYS:JD_end;


if verbose
   fprintf('Running Main loop %3.0f times\n', length(time_opt_vec));   
   fprintf('_______________________________________________________________________ \n');
   fprintf('| Iter.|         Date         |Delta-V|Duration|     Arrival Date     |\n')
   fprintf('_______________________________________________________________________ \n');
end
%% Main Optimization Loop 



for loop_id = 1:length(time_opt_vec)
    
    % Obtain Current time
    JD_loop = time_opt_vec(loop_id);
    
    % obtain the position of earth
    [pos_launch_eci, vel_launch_eci, ~] = ECEFtoECI(JD_loop,launch_pos_ecef',0,[0,0,0]');

    launch_state = [pos_launch_eci(:); vel_launch_eci(:)];
    
    % Obtain the position and velocity or primary and secondary
    [pos, vel] = planetEphemeris(JD_loop, primary_body, secondary_body);
    
    % Generate the ber_C_bei mat
    bei_C_ber = gen_beiCber_fun(pos,vel);
    % Need to loop with a function
    
    % fun_loop_opt(Cmat, convert, pos_earth,vel_earth, manifolds)
    [delta_vtot_inloop,trans_duration_inloop, mani_duration_inloop,...
        out_info_inloop] = manifold_transfer_lookup_fun(launch_state, ...
        manifolds, bei_C_ber, convert);
    
    % --- TODO --- Do moon constraint check for SEL2 HERE ----
    
    if delta_vtot_inloop <= delta_v_best
        delta_v_best = delta_vtot_inloop;
        t_launch_best = JD_loop;
        t_duration_best = trans_duration_inloop + mani_duration_inloop;
        out_info_best = out_info_inloop;
        out_info_best.launch_state = launch_state;
        
    end
    
    if 1
        
    total_dur_inloop = trans_duration_inloop + mani_duration_inloop;
    fprintf('| %4.0f | %s | %5.2f | %6.2f | %s |\n',...
        loop_id, datetime(JD_loop,'convertfrom', 'juliandate'),...
        delta_vtot_inloop/KILOMETERS, (total_dur_inloop)/DAYS,...
        datetime(JD_loop +total_dur_inloop/DAYS,'convertfrom', 'juliandate'))
    end
end

%% Check Moon Constraint (Not hitting the moon!)

%% Prepare and Return Output Output
delta_v = delta_v_best + atm_drag;

% convert time to time vector
t_launch = datevec(datetime(t_launch_best,'convertfrom', 'juliandate'));

t_duration = t_duration_best;

% out_struct should contain all information
out_struct = out_info_best;
out_struct.delta_v = delta_v;
out_struct.t_launch = t_launch;
out_struct.t_duration = t_duration;

out_struct.input_struct = input_struct;
out_struct.convert = convert;
out_struct.manifolds = manifolds;

if verbose
   fprintf('_______________________________________________________________________ \n');
   fprintf('Completed! \n\n');   
   
   fprintf('Best Launch Window: %s \n', datetime(t_launch_best,'convertfrom', 'juliandate'));
   fprintf('Total Transport Delta-V: %5.2f (km/s) \n', delta_v/KILOMETERS);
   fprintf('Total Transfer Duration: %6.2f (Days) \n', t_duration/DAYS);
   fprintf('Manifold Arrival Time: %s \n',  datetime(t_launch_best +t_duration/DAYS,'convertfrom', 'juliandate'));
end

end