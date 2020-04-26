
% performs simulation of one day attempting to get to a certain point on
% orbit.

% Assumptions
% - Launch is from cape canaveral
% - Azimuth range is limited
% - Assume ECEF and ECI is aligned at time = 0

close all;

%%  Universal values

font_size = 22;

flag_plot = 1;
save_flag = 0;

% universal values
mu = 3.986004418*10^(14)*METERS^3/SECONDS^2;
earth_radius = 6371*KILOMETERS;
v_earth = 491*METERS/SECONDS;
ecliptic_inc = -23.5*DEGREES;

a_SEL2 = 1.4959787e8 * KILOMETERS;

% solver options
opts = odeset('RelTol',1e-10,'AbsTol',1e-10);

%% Load data
load('results_2_18.mat')

% target_pos_0 = xmanifold.xnumvec_eci_m(1:3, 1);
% target_vel_0 = xmanifold.xnumvec_eci_m(4:6, 1);


%% Rotation Matrices for frame
% rotation matrix helper
cei_R_eciday = [1 0 0;
    0 cos(ecliptic_inc) -sin(ecliptic_inc);
    0 sin(ecliptic_inc) cos(ecliptic_inc)];

% sets the intiial offset for time
day_0 = 0; hour_0 = 0;
n_day = (2*pi/(23.93*HOURS));
n_year = (2*pi/(365*DAYS));
theta_hour = (day_0 + hour_0); % defines the initial offset

eciday_R_ecef = @(t) [cos(n_day*(theta_hour + t)) -sin(n_day*(theta_hour + t)) 0;
    sin(n_day*(theta_hour + t)) cos(n_day*(theta_hour + t)) 0;
    0 0 1];
eci_R_cei = @(t) [cos(n_year*(theta_hour + t)) -sin(n_year*(theta_hour + t)) 0;
    sin(n_year*(theta_hour + t)) cos(n_year*(theta_hour + t)) 0;
    0 0 1];
%% Location of cape canaveral in ECEF frame

% Longitude/Lattitude
lat = 28*DEGREES;
long = 20*DEGREES;

% Convert position of latitude longitude to ECEF frame
pos_ecef = [earth_radius*cos(lat)*cos(long); earth_radius*cos(lat)*sin(long);
    earth_radius*sin(lat)];


%% evolution of time
timevec = linspace(0*DAYS, 365*DAYS, 365*24); % time is an hours

day = 0*DAYS + 0*365.25/2*DAYS; % days
hour = 0*HOURS; % hours

results.timevec = timevec;
lambert_results.launch_pos = zeros(3, length(timevec));
lambert_results.launch_vel = zeros(3, length(timevec));
lambert_results.deltav1 = zeros(3, length(timevec));
lambert_results.deltav2 = zeros(3, length(timevec));
lambert_results.total_deltav = zeros(1, length(timevec));
lambert_results.collision = zeros(1, length(timevec));
lambert_results.target_pos = zeros(3, length(timevec));
lambert_results.target_vel = zeros(3, length(timevec));
lambert_results.solind_min = zeros(1,length(timevec));
lambert_results.t_manifold = zeros(1,length(timevec));
lambert_results.t_transfer = zeros(1,length(timevec));

figure;
hold on
plot_earth(earth_radius, 0);
view(3)
grid on
axis equal
xlabel('$x$ [$R_e$]' , 'Interpreter', 'Latex', 'FontSize', font_size)
ylabel('$y$ [$R_e$]' , 'Interpreter', 'Latex', 'FontSize', font_size)
zlabel('$z$ [$R_e$]' , 'Interpreter', 'Latex', 'FontSize', font_size)

for t_ind = 1:length(timevec)
    
    if mod(t_ind, 20) == 0      
        t_ind/(365*24)
    end
    curr_time = timevec(t_ind);
    
    [de, dn] = earth_vel_fun(pos_ecef);
    
    % Obtain the position based on rotation of the day
    pos_launch_eciday = eciday_R_ecef(curr_time)*pos_ecef;
    vel_launch_eciday = eciday_R_ecef(curr_time)*v_earth*de;
    
    % Obtain the correction with the inclination
    pos_launch_cei = cei_R_eciday*pos_launch_eciday;
    vel_launch_cei = cei_R_eciday*vel_launch_eciday;
    
    % This is not an ECI frame this is an CEI frame (it is misnamed)
    pos_launch_eci = eci_R_cei(curr_time)*pos_launch_cei;
    vel_launch_eci = eci_R_cei(curr_time)*vel_launch_cei;
    
    % ---------- Lamberts problem ----------
    % begin min search
    total_deltav_min = Inf;
    solind_min = NaN;
    v1_lamvec_min = NaN;
    v2_lamvec_min = NaN;
    
    target_pos_min = zeros(3,1);
    target_vel_min = zeros(3,1);
    for ii = 2:20
        
        orb_ind = results.orb_sort_ind(ii);
       
        
        target_pos_0 = results.sol{orb_ind}.xnumvec_cei(1:3,1)*a_SEL2;
        target_vel_0 = results.sol{orb_ind}.xnumvec_cei(4:6,1)*a_SEL2/YEARS*2*pi;
        
        target_pos = eci_R_cei(curr_time)*target_pos_0;
        target_vel = eci_R_cei(curr_time)*target_vel_0;
        
        x1_state = [pos_launch_eci; vel_launch_eci];
        x2_state = [target_pos; target_vel];
        % minimum a lambert problem solution
        [v1_lamvec,v2_lamvec, T_lamb] = lambert_mina(mu, x1_state, x2_state);
        
        % general lambert problem solution
        % [vi, vf] = glambert(mu, x1_state, x2_state, 5*T_lamb, -1);
        
        xinit_lamb = [pos_launch_eci; v1_lamvec];
        
        [t_trans, xvec_trans] = propagate_lambert_fun(mu, T_lamb, xinit_lamb, opts);
        
        % check for collision
        col = check_collision_fun(xvec_trans(:,1:3), 0.95*earth_radius);
        
        % ------------- lambert_results ----------
        
        total_deltav = norm(v1_lamvec-vel_launch_eci) + norm(v2_lamvec - target_vel);
        
        if total_deltav <= total_deltav_min && ~col
            T_lamb_min = T_lamb/DAYS;
            total_deltav_min = total_deltav;
            solind_min = orb_ind;
            t_min_manifold = -results.sol{solind_min}.t(results.minind(solind_min))/(2*pi)*YEARS/DAYS;
%             ii
            v1_lamvec_min = v1_lamvec;
            v2_lamvec_min = v2_lamvec;
            target_pos_min = target_pos;
            target_vel_min = target_vel;
            
            if flag_plot && ~col
                % ----------- plotting ----------------
                p1 = plot3(pos_launch_eci(1),pos_launch_eci(2),pos_launch_eci(3), 'ko');
                quiver3(pos_launch_eci(1),pos_launch_eci(2),pos_launch_eci(3),1e3*vel_launch_eci(1), 1e3*vel_launch_eci(2), 1e3*vel_launch_eci(3),'k-', 'linewidth', 2)
                p2 = plot3(xvec_trans(:,1), xvec_trans(:,2), xvec_trans(:,3), 'k--');
                p3 = plot3(target_pos(1),target_pos(2),target_pos(3), 'k*');
                p4 = quiver3(target_pos(1),target_pos(2),target_pos(3),3e2*target_vel(1), 3e2*target_vel(2), 3e2*target_vel(3),'k-', 'linewidth', 2);
                drawnow limitrate
                view([-54 30])
                pause(0.5)
            end
            
        end
    end
    
    taget_pos = target_pos_min;
    taget_vel = target_vel_min;
    
    if total_deltav_min == Inf
        col = 1;
    end
    
    % save values
    lambert_results.launch_pos(:,t_ind) = pos_launch_eci;
    lambert_results.launch_vel(:,t_ind) = vel_launch_eci;
    lambert_results.deltav1(:,t_ind) = v1_lamvec_min;
    lambert_results.deltav2(:,t_ind) = v2_lamvec_min;
    lambert_results.solind_min(t_ind) = solind_min;
    lambert_results.target_pos(:,t_ind) = target_pos_min;
    lambert_results.target_vel(:,t_ind) = target_vel_min;
    lambert_results.t_manifold(t_ind) = t_min_manifold;
    lambert_results.t_transfer(t_ind) = T_lamb_min;

    
    if col
        lambert_results.total_deltav(:,t_ind) = NaN;
    else % if no collision all results exist
        lambert_results.total_deltav(:,t_ind) = total_deltav_min;
        
    end
    
    lambert_results.collision(t_ind) = col;
    
    if flag_plot && ~col
        % ----------- plotting ----------------
        p1 = plot3(pos_launch_eci(1),pos_launch_eci(2),pos_launch_eci(3), 'ko');
        quiver3(pos_launch_eci(1),pos_launch_eci(2),pos_launch_eci(3),1e3*vel_launch_eci(1), 1e3*vel_launch_eci(2), 1e3*vel_launch_eci(3),'k-', 'linewidth', 2)
        p2 = plot3(xvec_trans(:,1), xvec_trans(:,2), xvec_trans(:,3), 'k--');
        p3 = plot3(target_pos(1),target_pos(2),target_pos(3), 'k*');
        p4 = quiver3(target_pos(1),target_pos(2),target_pos(3),3e2*target_vel(1), 3e2*target_vel(2), 3e2*target_vel(3),'k-', 'linewidth', 2);
        drawnow limitrate
        view([-54 30])
        pause(0.5)
    end
    
    % save just in case
    if mod(t_ind, 1000) == 0 && save_flag
        fprintf('Saving info ...')
        save('lambert_results.mat', 'lambert_results')
        fprintf('done! \n')
    end
    %     delete(p)
end


function [de, dn] = earth_vel_fun(pos_launch_vec)

% north pole in ECEF
n_ecef = [0 0 1]'; % nort h pole

a_vec = pos_launch_vec/norm(pos_launch_vec);

de = cross(n_ecef, a_vec); % tangential due east vector

% due east vector
de = de/norm(de);

% due north vector  from launch
dn = cross(a_vec, de);

dn = dn/norm(dn);
end

function [t_trans, xvec_trans] = propagate_lambert_fun(mu, T_lamb, xinit_lamb, opts)

[t_trans,xvec_trans] = ode45(@(t,x) odekep(t,x,0,mu,1), ...
    [0 T_lamb], xinit_lamb, opts);

end

function col = check_collision_fun(xvec, earth_radius)

dist_to_center = sum(xvec.*xvec,2);

min_dist = min(dist_to_center - earth_radius^2);

if min_dist < 0
    
    col = 1;
else
    col = 0;
    
end

end