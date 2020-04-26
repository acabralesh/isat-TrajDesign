% hohman transfer script

% test the hohman transfer

clc; 
% close all; 

addpath('Utils')
addpath('../Units')

warning off

%% Global Parameters
font_size = 18;

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

% Velocity of earth

% non dimensionalized system rotation
omega_rot = [0, 0, 1]';

% non dimensionalized velocity of earth
v_earth_inert = cross(omega_rot, [1-mu, 0, 0]');

% ODE sets
opts = odeset('RelTol',1e-14,'AbsTol',1e-14);
%% Find Location of L1 L2 based on roots of polinomial 
L1 = roots([1,-(3-mu),(3-2*mu),-mu,2*mu,-mu]);
L2 = roots([1,(3-mu),(3-2*mu),-mu,-2*mu,-mu]);

% x-position of Lagrange points in rotating frame
% Note xL1 is non dimensional
x_l1 = L1(imag(L1)==0 & real(L1)>0);
x_l2 = L2(imag(L2)==0 & real(L2)>0);

% desired point:
x_l = x_l2;

%% LEO orbit

leo_alt = 400*KILOMETERS; %[m]

leo_r_mag = 6371*KILOMETERS + leo_alt;  %[m]

leo_vel_mag = sqrt(G*M_earth/leo_r_mag); %[m/s]

% non dimensionalizing units
leo_r_mag = leo_r_mag/a;

leo_vel_mag = leo_vel_mag/a*YEARS/(2*pi);

leo_r_earth = [-leo_r_mag; 0; 0];

leo_v_earth_inert = [0; -leo_vel_mag; 0];

leo_t_period = 2*pi*sqrt(leo_r_mag^3/mu);

% generate the initial condition in baricenter rotating frame
[r_bari_inert, v_orb_rot_sol] = convert_bari_rot(1-mu, ...
    v_earth_inert, omega_rot, leo_r_earth, leo_v_earth_inert);
x_leo_init = [r_bari_inert;v_orb_rot_sol];

%% L2 orbit
l2_r_mag = (x_l2)*a; % [m]

l2_vel_mag = sqrt(G*M_earth/l2_r_mag); %[m/s]

% non dimensionalizing units
l2_r_mag = l2_r_mag/a;

l2_vel_mag = l2_vel_mag/a*YEARS/(2*pi);

l2_r_earth = [l2_r_mag; 0; 0];

l2_v_earth_inert = [0; l2_vel_mag; 0];

l2_t_period = 2*pi*sqrt(l2_r_mag^3/mu);

[r_bari_inert, v_orb_rot_sol] = convert_bari_rot(1-mu, ...
    v_earth_inert, omega_rot, l2_r_earth, l2_v_earth_inert);
x_l2_init = [r_bari_inert;v_orb_rot_sol];


%% Computing hohman transfer

r1 = leo_r_mag;
r2 = l2_r_mag;

deltav1 = sqrt(mu/r1)*(sqrt(2*r2/(r1+r2)) - 1);

% minimum value should be the on the same order as a circularize + a lil
% bit more, difference is that although you dont have to circularize it
% then it will be better to launch as close as possible 

trans_t_period = 2*pi*sqrt(((r1 +r2)/2)^3/mu);

% angle of shooting wrt to -x axis
ang = 0*pi/8 ; % correct angle to hit L2 = pi/6.65;
R = [cos(ang), -sin(ang) 0; sin(ang) cos(ang) 0; 0 0 1];

trans_r_earth = R*leo_r_earth;

trans_v_earth = R*(leo_v_earth_inert + [0; -deltav1; 0]);

[r_bari_inert, v_orb_rot_sol] = convert_bari_rot(1-mu, ...
    v_earth_inert, omega_rot, trans_r_earth, trans_v_earth);
x_trans_init = [r_bari_inert;v_orb_rot_sol];

trans_t_period = 29.6*DAYS/YEARS*2*pi;
%% Compute LEO L2 orbit

[t,leo_xnumvec] = ode45(@(t,x) eom_3body_rot_fun(t,x,mu), ...
    [0 leo_t_period], x_leo_init, opts);

x_l2_init(5) = 0;
[t,l2_xnumvec] = ode45(@(t,x) eom_3body_rot_fun(t,x,mu), ...
    [0 l2_t_period], x_l2_init, opts);

% transfer
[t,leotrans_xnumvec] = ode45(@(t,x) eom_3body_rot_fun(t,x,mu), ...
    [0 trans_t_period], x_trans_init, opts);

figure;
hold on
% plot Earth
% --- earth ------
r_radius = 6371*KILOMETERS/a;
x_offset = (1-mu);
plot_earth(r_radius, x_offset)

% LEO
plot3(leo_xnumvec(:,1), leo_xnumvec(:,2), leo_xnumvec(:,3), '--k', ...
    'linewidth', 1);

% L2
plot3(l2_xnumvec(:,1), l2_xnumvec(:,2), l2_xnumvec(:,3), '-*k', ...
    'linewidth', 1);

% trans
plot3(leotrans_xnumvec(:,1), leotrans_xnumvec(:,2), leotrans_xnumvec(:,3), '-k', ...
    'linewidth', 1);

% axis equal
xlabel('$x$ [Au]' , 'Interpreter', 'Latex', 'FontSize', font_size)
ylabel('$y$ [Au]' , 'Interpreter', 'Latex', 'FontSize', font_size)
zlabel('$z$ [Au]' , 'Interpreter', 'Latex', 'FontSize', font_size)
set(gca,'fontsize',font_size-4)
% ylim([-1e-3 1e-3])
grid on


%% Converting to rotating frame

function [r_bari_inert, v_orb_rot_sol] = convert_bari_rot(dist_2_earth, ...
    v_earth_inert, omega_rot, r_earth, v_orb_inert_earth)

r_bari_inert = [dist_2_earth; 0; 0] + r_earth;

v_orb_inert_sol = v_orb_inert_earth + v_earth_inert;

v_orb_rot_sol = v_orb_inert_sol - cross(omega_rot, r_bari_inert);

end
%% Plotting Earth
function plot_earth(r_radius, xoffset)
[xearth, yearth, zearth] = sphere;
hSurface=surf(xoffset + xearth*r_radius,yearth*r_radius,zearth*r_radius);
set(hSurface,'FaceColor',[0 0 1], ...
      'FaceAlpha',0.5,'FaceLighting','gouraud')
end

