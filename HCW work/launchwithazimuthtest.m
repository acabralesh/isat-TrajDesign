% Introduction of ECI frame into the problem
close all;
addpath('Utils')
addpath('../Units')
font_size = 22;
% universal values

mu = 3.986004418*10^(14)*METERS^3/SECONDS^2;

earth_radius = 6371*KILOMETERS;

v_earth = 0*491*METERS/SECONDS;

opts = odeset('RelTol',1e-10,'AbsTol',1e-10);


% ECEF frame
lat = 28*DEGREES;
long = 20*DEGREES;

bearing = 90*DEGREES;

% Convert position of latitude longitude to ECEF frame


pos_ecef = [earth_radius*cos(lat)*cos(long); earth_radius*cos(lat)*sin(long);
    earth_radius*sin(lat)];

n_ecef = [0 0 1]'; % nort h pole

a_vec = pos_ecef/norm(pos_ecef);

de = cross(n_ecef, a_vec); % tangential due east vector

% due east vector
de = de/norm(de);

% due north vector  from launch
dn = cross(a_vec, de);

dn = dn/norm(dn);

% velocity vector given the direction
d = dn*cos(bearing) + de*sin(bearing);

% Obtain orbit requirements
r_orb = earth_radius + 500*KILOMETERS;

dv_launch = sqrt(mu/earth_radius)*(sqrt(2*r_orb/(r_orb + earth_radius)));


xinit = [pos_ecef; dv_launch*d];

T = 2*pi*sqrt(r_orb^3/mu);

t = [0 2*T];


[t,xvec_sc1] = ode45(@(t,x) odekep(t,x,0,mu,1), ...
        t, xinit, opts);


t = [0 8*T];
t = [0 2.88*HOURS];
figure;
for ii = 1:100:5000

load('xmanifold.mat')    
t = xmanifold.t;
xinit_mani = xmanifold.xnumvec_eci_m(1:6,ii);

[t,xvec_mani] = ode45(@(t,x) odekep(t,x,0,mu,1), ...
        t, xinit_mani, opts);



plot_earth(1, 0);
hold on
plot3(xvec_sc1(:,1)/earth_radius, xvec_sc1(:,2)/earth_radius, xvec_sc1(:,3)/earth_radius,'k-', 'linewidth', 1)
quiver3(a_vec(1), a_vec(2), a_vec(3), 0.4*de(1), 0.4*de(2), 0.4*de(3),'k-', 'linewidth', 1)
quiver3(a_vec(1), a_vec(2), a_vec(3), 0.4*dn(1), 0.4*dn(2), 0.4*dn(3),'k-', 'linewidth', 1)
p = plot3(xvec_mani(:,1)/earth_radius,xvec_mani(:,2)/earth_radius,xvec_mani(:,3)/earth_radius) 
plot3(xmanifold.xnumvec_eci_m(1,:)/earth_radius,xmanifold.xnumvec_eci_m(2,:)/earth_radius,xmanifold.xnumvec_eci_m(3,:)/earth_radius) 
d = plot3(xinit_mani(1)/earth_radius,xinit_mani(2)/earth_radius,xinit_mani(3)/earth_radius, '*k') 

axis square
xlabel('$x$ [$R_e$]' , 'Interpreter', 'Latex', 'FontSize', font_size)
ylabel('$y$ [$R_e$]' , 'Interpreter', 'Latex', 'FontSize', font_size)
zlabel('$z$ [$R_e$]' , 'Interpreter', 'Latex', 'FontSize', font_size)
set(gca,'fontsize',font_size-4)
grid on
axis equal
% xlim([-5 5])
% ylim([-5 5])
% zlim([-5 5])
pause(0.1)

delete(p)
delete(d)

end
[a, eccentricity, inclination, longnode, argperi, meananom] = rec_to_kepler(mu, xinit(1:3), xinit(4:6));

inclination/DEGREES

acos(sin(bearing)*cos(lat))/DEGREES