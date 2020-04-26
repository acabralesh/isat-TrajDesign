% Testing Script
close all
addpath('../Units')
addpath('Utils')
font_size = 18;

mu = 3.986004418*10^(14)*METERS^3/SECONDS^2;

earth_radius = 6351*KILOMETERS;

SIM.UNI.mu = mu;
SC.Phy.mass = 100*KILOGRAMS;

force = 0*ones(3,1);

tspan = [0:1*MINUTES:2*HOURS];

r_orb_mag = earth_radius + 1000*KILOMETERS;

v_orb_mag = sqrt(mu/r_orb_mag);

x0 =[r_orb_mag, 0, 0, 0, v_orb_mag, 0]';
n = v_orb_mag/r_orb_mag;
 
xr = [5*METERS, -10*METERS, 2, -5*n,-2*5*n,2*n]';

% xr = [5*METERS, 0, 0, 0,-2*n*5,0]';

[SC2_pos_ECI,SC2_vel_ECI] = lvlh_to_eci(x0(1:3), x0(4:6),xr(1:3),xr(4:6));


x1 =[SC2_pos_ECI; SC2_vel_ECI];

opts = odeset('RelTol',1e-10,'AbsTol',1e-10);

[t,x0] = ode45(@(t,x) odekep(t,x,force,SIM,SC), tspan, x0, opts);

[t,x1] = ode45(@(t,x) odekep(t,x,force,SIM,SC), tspan, x1, opts);

% Orbit
figure;
sphere;
hold on;
plot3(x0(:,1)/earth_radius, x0(:,2)/earth_radius,x0(:,3)/earth_radius, 'k--', 'linewidth', 2)
plot3(x1(:,1)/earth_radius, x1(:,2)/earth_radius,x1(:,3)/earth_radius, 'b--', 'linewidth', 2)
xlabel('$x$ $[R_e]$', 'Interpreter', 'Latex', 'FontSize', font_size);
ylabel('$y$ $[R_e]$', 'Interpreter', 'Latex','FontSize', font_size);
zlabel('$z$ $[R_e]$', 'Interpreter', 'Latex','FontSize', font_size);
set(gca,'fontsize',font_size-2)
axis equal
grid on

% Evolution of state over time
figure
plot(t/HOURS, x0/earth_radius)
xlabel('$t$ [Hours]', 'Interpreter', 'Latex','FontSize', font_size);
xlabel('[$R_e$]', 'Interpreter', 'Latex','FontSize', font_size);
title('Evolution of states', 'Interpreter', 'Latex', 'FontSize', font_size)
leg = legend({'$x$','$y$', '$z$', '$v_x$', '$v_y$', '$v_z$'}, 'Interpreter', 'Latex'); 
leg.FontSize = font_size;
set(gca,'fontsize',font_size-2)
grid on

% convert to LVLH frame
x1_lvlh = zeros(size(x1));

for ii = 1:length(x1)
   [pos1_temp, vel1_temp, a1_temp, ~] = eci_to_lvlh(mu, ...
    x0(ii,1:3)', x0(ii,4:end)', x1(ii,1:3)', x1(ii,4:end)');

x1_lvlh(ii,:) = [pos1_temp', vel1_temp'];
    
end

figure;
subplot(2,2,1)
plot(x1_lvlh(:,1)/METERS, x1_lvlh(:,2)/METERS, 'k--', 'linewidth', 1)
xlabel('$x$ [m]', 'Interpreter', 'Latex','FontSize', font_size);
ylabel('$y$ [m]', 'Interpreter', 'Latex','FontSize', font_size);
set(gca,'fontsize',font_size-2)
grid on
subplot(2,2,2)
plot(x1_lvlh(:,2)/METERS, x1_lvlh(:,3)/METERS, 'k--', 'linewidth', 1)
xlabel('$y$ [m]', 'Interpreter', 'Latex','FontSize', font_size);
ylabel('$z$ [m]', 'Interpreter', 'Latex','FontSize', font_size);
set(gca,'fontsize',font_size-2)
grid on
subplot(2,2,3)
plot(x1_lvlh(:,1)/METERS, x1_lvlh(:,3)/METERS, 'k--', 'linewidth', 1)
xlabel('$x$ [m]', 'Interpreter', 'Latex','FontSize', font_size);
ylabel('$z$ [m]', 'Interpreter', 'Latex','FontSize', font_size);
set(gca,'fontsize',font_size-2)
grid on
subplot(2,2,4)
plot3(x1_lvlh(:,1)/METERS,x1_lvlh(:,2)/METERS, x1_lvlh(:,3)/METERS, 'k--', 'linewidth', 1)
xlabel('$x$ [m]', 'Interpreter', 'Latex','FontSize', font_size);
ylabel('$y$ [m]', 'Interpreter', 'Latex','FontSize', font_size);
zlabel('$z$ [m]', 'Interpreter', 'Latex','FontSize', font_size);
set(gca,'fontsize',font_size-2)
grid on


function [xdot] = odekep(t,x,force, SIM, SC)

xdotkep = kepler_accel(SIM, SC, x, force);

xdot = [x(4:6); xdotkep];

end


