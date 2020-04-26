% unit test for the code
close all

%% Unit test that CEI is working

% Longitude/Lattitude
lat = 20*DEGREES;
long = 20*DEGREES;

% Convert position of latitude longitude to ECEF frame
pos_ecef = [earth_radius*cos(lat)*cos(long); earth_radius*cos(lat)*sin(long);
    earth_radius*sin(lat)];

% eciday_R_ecef used to spin during day
% eci_R_cei used to spin during time

time_1 = 0; % time zero, should be aligned with x axis
time_2 = 91*DAYS; % 12 hours later

time_3 = 180*DAYS; % 12 hours later

time_4 = 1*YEARS; % 12 hours later

timevec = [time_1, time_2, time_3,time_4];


r_eci = zeros(3,length(timevec));

figure;
hold on
plot_earth(earth_radius, 0);
view(3)
grid on
axis equal
xlabel('$x$ [$R_e$]' , 'Interpreter', 'Latex', 'FontSize', font_size)
ylabel('$y$ [$R_e$]' , 'Interpreter', 'Latex', 'FontSize', font_size)
zlabel('$z$ [$R_e$]' , 'Interpreter', 'Latex', 'FontSize', font_size)

for ii = 1:length(timevec)
    curr_time = timevec(ii);
    
% Obtain the position based on rotation of the day
pos_launch_eciday = eciday_R_ecef(curr_time)*pos_ecef;
% vel_launch_eciday = eciday_R_ecef(curr_time)*v_earth*de;
    
% Obtain the correction with the inclination
pos_launch_cei = cei_R_eciday*pos_launch_eciday;
% vel_launch_cei = cei_R_eciday*vel_launch_eciday;
    
% This is not an ECI frame this is an CEI frame (it is misnamed)
pos_launch_eci = eci_R_cei(curr_time)*pos_launch_cei;
% vel_launch_eci = eci_R_cei(curr_time)*vel_launch_cei;

r_eci(:,ii) = pos_launch_eci;

p1 = plot3(pos_launch_eci(1),pos_launch_eci(2),pos_launch_eci(3), 'ko');

end

%% Convert orbits to ECI frame

n_manifolds = 12;

for ii = 2:n_manifolds

    
end

