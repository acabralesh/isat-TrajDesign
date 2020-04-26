% test script
% close all

addpath('Utils')
addpath('../Units')
%% Global Parameters
font_size = 22;

% my values
G = 6.67408e-11 * METERS^3 * KILOGRAMS^-1 * SECONDS^-2;
M_sun = 1.9884754e30 * KILOGRAMS;
M_earth = 5.972358e24 * KILOGRAMS;
M_moon = 7.34767309e22 * KILOGRAMS;

% Sun - Earth
a = 1.4959787e8 * KILOMETERS;
% Earth - Moon
% a = 384400 * KILOMETERS;

% Sun - Earth
mu = M_earth/(M_sun + M_earth);

opts = odeset('RelTol',1e-14,'AbsTol',1e-14);


load('results.mat')

% modified
figure
for dx = 1:365
%% TIME

omega_n_ber = [0; 0; 1];

n = omega_n_ber(3); 

ecliptic_inc = -23.5*DEGREES;

day = dx*DAYS + 0*365.25/2*DAYS; % days

hour = 0*HOURS; % hours

theta_year = n* (day + hour)* (2*pi/YEARS); % converts days to non dimensional

theta_hour = n*(2*pi/DAYS)*(day + hour);

bei_R_ber = @(t) [cos(theta_year + t) -sin(theta_year + t) 0; 
   sin(theta_year + t) cos(theta_year + t) 0; 
   0 0 1];

eci_R_ecef = @(t) [cos(theta_hour + t) -sin(theta_hour + t) 0; 
   sin(theta_hour + t) cos(theta_hour + t) 0; 
   0 0 1];


cei_R_eci = [1 0 0;
            0 cos(ecliptic_inc) -sin(ecliptic_inc);
            0 sin(ecliptic_inc) cos(ecliptic_inc)];
              
%% 


testid = results.orb_sort_ind(4);

xnumvec_ber = results.sol{testid}.xnumvec';

t_min_ind = results.minind(testid);

% obtain the minimum time
tmin = results.sol{testid}.t(t_min_ind);

% obtain the manifold from start to finish
xnumvec_ber = xnumvec_ber(:, t_min_ind:-1:1);
time = results.sol{testid}.t(t_min_ind:-1:1) - tmin;

% ------ % test 
[sc_pos_bei,sc_vel_bei] = ber_to_bei(bei_R_ber(0), omega_n_ber,xnumvec_ber(1:3,:), xnumvec_ber(4:6,:));

cei_R_bei = eye(3);
r_earth_bei = bei_R_ber(0)*[1-mu;0;0];

[sc_pos_cei,sc_vel_cei] = bei_to_cei(cei_R_bei, omega_n_ber,r_earth_bei, sc_pos_bei, sc_vel_bei);

% figure
% plot3(sc_pos_cei(1,:)*a/earthRadius,sc_pos_cei(2,:)*a/earthRadius,...
%     sc_pos_cei(3,:)*a/earthRadius)
% axis equal
% grid on
% --------
xnumvec_bei = zeros(size(xnumvec_ber));
xnumvec_cei = zeros(size(xnumvec_ber));
xnumvec_eci = zeros(size(xnumvec_ber));


for ii = 1:length(xnumvec_ber)
    
    % obtain time
    theta = time(ii);
    
    % Convert BER to BEI
    [sc_pos_bei,sc_vel_bei] = ber_to_bei(bei_R_ber(theta), omega_n_ber,xnumvec_ber(1:3,ii), xnumvec_ber(4:6,ii));
    
    xnumvec_bei(:, ii) = [sc_pos_bei; sc_vel_bei];
    
    % Convert BER to CEI
    r_earth_bei = bei_R_ber(theta)*[1-mu;0;0];
    
    [sc_pos_cei,sc_vel_cei] = bei_to_cei(cei_R_bei, omega_n_ber,r_earth_bei, sc_pos_bei, sc_vel_bei);

    xnumvec_cei(:, ii) = [sc_pos_cei;sc_vel_cei];
    
    sc_pos_eci = cei_R_eci'*sc_pos_cei;
    sc_vel_eci = cei_R_eci'*sc_vel_cei;
    
    xnumvec_eci(:, ii) = [sc_pos_eci;sc_vel_eci];
    
end

% figure(1)
% plot3(xnumvec_bei(1,:),xnumvec_bei(2,:), xnumvec_bei(3,:))
% axis equal
% grid on

% figure
hold on
plot_earth(1, 0)
p = plot3(xnumvec_cei(1,:)*a/earthRadius,xnumvec_cei(2,:)*a/earthRadius,...
    xnumvec_cei(3,:)*a/earthRadius)
d = plot3(xnumvec_eci(1,:)*a/earthRadius,xnumvec_eci(2,:)*a/earthRadius,...
    xnumvec_eci(3,:)*a/earthRadius, 'k')
drawnow 

axis equal
xlim([-10 10])
ylim([-10 10])
zlim([-10 10])
% view([42 17])
grid on

xlabel('$x$ [$R_e$]' , 'Interpreter', 'Latex', 'FontSize', font_size)
ylabel('$y$ [$R_e$]' , 'Interpreter', 'Latex', 'FontSize', font_size)
zlabel('$z$ [$R_e$]' , 'Interpreter', 'Latex', 'FontSize', font_size)


% pause(0.1)
delete(p)
delete(d)
end
