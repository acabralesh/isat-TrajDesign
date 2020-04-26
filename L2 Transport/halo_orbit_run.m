%% Solving for HALO orbits

clc; 
close all; 

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

% mu = M_moon/(M_moon + M_earth);
% my values
%% Desired Halo Parameters

Az = 110000 * KILOMETERS; % kilometers

phi = 0;

n = 1; % 1 for class 1, 3 for class 2

libind = 2;
% 
tol = 1e-12;

%% Finding Libration points

[x_l1,x_l2] = find_librationpts(mu);

% desired point:
x_l = [x_l1,x_l2];

x_l = x_l(libind);

%% Setting up Richardson's Solution

gamma_l = x_l;

Az = Az/(gamma_l*a);
rE = gamma_l*a; % distance between primary to libration point

% Obtain the initial conditions for 3rd order solution
[x0vec] = halo_3bodyrot_richardson_init(mu, x_l, libind, Az, phi, n);

t = 0:0.01:2*pi;

% obtain the vectorize solution for one orbit of libration point
[x, y, z, xd, yd, zd] = halo_3bodyrot_richardson_vec(...
    t, mu, x_l, libind, Az, phi, n);

xlin = x*rE;
ylin = y*rE;
zlin = z*rE;

figure
plot3(x, y, z,'k-')
xlabel('$x$ [$\gamma_L$]' , 'Interpreter', 'Latex', 'FontSize', font_size)
ylabel('$y$ [$\gamma_L$]' , 'Interpreter', 'Latex', 'FontSize', font_size)
zlabel('$z$ [$\gamma_L$]' , 'Interpreter', 'Latex', 'FontSize', font_size)
set(gca,'fontsize',font_size-4)
grid on
axis equal
view([-41 21])

%% Obtaining the solution through differential correctors

% set up the correct initial condition in the same units
libchange = [-1, 1];

% x0vec is in frame relative to L2, need to change it to baricenter frame
x0_ = (x0vec*rE + [1 - mu + libchange(libind)*gamma_l, 0, 0, 0, 0, 0]'*a)/a;

% perform the differential corrector scheme
[x0_dc,t_per, num_iter] = diffcorr_halo_3bodyrot(mu, x0_, tol);

%% Test the differential corrector solution

opts = odeset('RelTol',1e-14,'AbsTol',1e-14);
% tspan = [0, 4.5*pi];
tspan = [0, t_per];

[t,xnumvec] = ode45(@(t,x) eom_3body_rot_fun(t,x,mu), tspan, x0_dc, opts);

% Halo in km
xa = (x*rE + (1-mu+gamma_l)*a)/a;
ya = y*rE/a;
za = z*rE/a;

xnum = xnumvec(:,1);
ynum = xnumvec(:,2);
znum = xnumvec(:,3);

figure
% Y Z Axis
subplot(2,2,1)
plot(ynum, znum,'k-')
hold on
plot(ya,za, 'b--')
xlabel('$y$ [Au]' , 'Interpreter', 'Latex', 'FontSize', font_size)
ylabel('$z$ [Au]' , 'Interpreter', 'Latex', 'FontSize', font_size)
set(gca,'fontsize',font_size-4)
grid on
axis square
% X Z axis
subplot(2,2,2)
plot(xnum, znum,'k-')
hold on
plot(xa, za,'b--')
xlabel('$x$ [Au]' , 'Interpreter', 'Latex', 'FontSize', font_size)
ylabel('$z$ [Au]' , 'Interpreter', 'Latex', 'FontSize', font_size)
set(gca,'fontsize',font_size-4)
grid on
axis square

% X Y axis
subplot(2,2,3)
plot(xnum, ynum,'k-')
hold on
plot(xa, ya,'b--')
xlabel('$x$ [Au]' , 'Interpreter', 'Latex', 'FontSize', font_size)
ylabel('$y$ [Au]' , 'Interpreter', 'Latex', 'FontSize', font_size)
set(gca,'fontsize',font_size-4)
grid on
axis square
% 3d
subplot(2,2,4)
plot3(xnum, ynum, znum,'k-')
hold on
plot3(xa, ya, za,'b--')
xlabel('$x$ [Au]' , 'Interpreter', 'Latex', 'FontSize', font_size)
ylabel('$y$ [Au]' , 'Interpreter', 'Latex', 'FontSize', font_size)
zlabel('$z$ [Au]' , 'Interpreter', 'Latex', 'FontSize', font_size)
set(gca,'fontsize',font_size-4)
grid on
axis square
view([-41 21])

geoalt = 42164*KILOMETERS/a;
moonalt = 384400*KILOMETERS/a;
geo_orb = [geoalt*cos(2*t), geoalt*sin(2*t), 0*cos(2*t)];
moon_orb = [moonalt*cos(3*t), moonalt*sin(3*t), 0*cos(2*t)];

figure; 
plot3(1, 0, 0, '*k')
hold on
plot3(1+ geo_orb(:,1), geo_orb(:,2), geo_orb(:,3), '--k')
plot3(1+ moon_orb(:,1), moon_orb(:,2), moon_orb(:,3), '--b')
plot3(xa, ya, za,'b--')
plot3(xnum, ynum, znum,'k-')
xlabel('$x$ [Au]' , 'Interpreter', 'Latex', 'FontSize', font_size)
ylabel('$y$ [Au]' , 'Interpreter', 'Latex', 'FontSize', font_size)
zlabel('$z$ [Au]' , 'Interpreter', 'Latex', 'FontSize', font_size)
set(gca,'fontsize',font_size-4)
grid on
view([-16, 5])
leg = legend({'Earth', 'GEO','Moon', 'linear halo', 'nonlinear halo'}, 'Interpreter', 'Latex');
leg.FontSize = font_size;
axis square

%% Manifolds

% time vector 

t_stm_span = 0:0.05:t_per;

% Initial Conditions for halo orbit
xstate_init = x0_dc;

STM_0 = eye(6);

% Separation along stable eigenvector;
d = 200*KILOMETERS/a;

t_mani_span = [0 -2.4*t_per];

sol = halo_find_stablemani(mu, t_stm_span,...
    xstate_init, STM_0, d, t_mani_span, opts);

%%
figure; 
plot3(1 - mu, 0, 0, 'ob', 'MarkerFaceColor', [0 0 1])
hold on
% plot3(1+ geo_orb(:,1), geo_orb(:,2), geo_orb(:,3), '--k')
% plot3(1+ moon_orb(:,1), moon_orb(:,2), moon_orb(:,3), '--b')
plot3(xnum, ynum, znum,'k-', 'linewidth', 2)
xlabel('$x$ [Au]' , 'Interpreter', 'Latex', 'FontSize', font_size)
ylabel('$y$ [Au]' , 'Interpreter', 'Latex', 'FontSize', font_size)
zlabel('$z$ [Au]' , 'Interpreter', 'Latex', 'FontSize', font_size)
set(gca,'fontsize',font_size-4)
grid on
view([-16, 5])

for ii = 1:length(sol)
   xnumvec = sol{ii}.xnumvec;
   
   plot3(xnumvec(:,1),xnumvec(:,2),xnumvec(:,3),'-','Color',[0.6 0.6 0.6], 'linewidth', 0.1)
% patchline(xnumvec(:,1),xnumvec(:,2),xnumvec(:,3),'linestyle','--','edgecolor','b','linewidth',2,'edgealpha',0.8);
end
patchline(1- mu+ geo_orb(:,1), geo_orb(:,2), geo_orb(:,3),'linestyle','-','edgecolor','k','linewidth',1,'edgealpha',0.8);
patchline(1 - mu + moon_orb(:,1), moon_orb(:,2), moon_orb(:,3),'linestyle','-','edgecolor',[0.3 0.3 0.3],'linewidth',1,'edgealpha',0.8);
% plot3(1 - mu + geo_orb(:,1), geo_orb(:,2), geo_orb(:,3), '-k', 'linewidth', 2)
% plot3(1 - mu + moon_orb(:,1), moon_orb(:,2), moon_orb(:,3), '-', 'linewidth', 2)
plot3(xnum, ynum, znum,'k-', 'linewidth', 2)

% leg = legend({'Earth', 'GEO','Moon', 'nonlinear halo'}, 'Interpreter', 'Latex');
% leg.FontSize = font_size;
axis equal
xlim([.992, 1.015]);
ylim([-7e-3 7e-3]);
zlim([-1e-3 1e-3]);

%% Finding the orbit closest to LEO

leoalt = (429+6371) * KILOMETERS/a;
leoalt = (327+6371) * KILOMETERS/a;
tz = 0:0.01:2*pi;
leo_orb = [leoalt*cos(tz'), leoalt*sin(tz'), 0*cos(tz')];

% orbit closest to LEO


% computes the minimum distance to earth and the corresponding indice
mindist = zeros(length(sol),1);
minind = zeros(length(sol),1);

for ii = 1:length(sol)
    xnumvec = sol{ii}.xnumvec;
    X = xnumvec(:,1) - (1 - mu);
    Y = xnumvec(:,2) - 0;
    Z = xnumvec(:,3) - 0;
    
    DIST = sqrt(X.^2 + Y.^2 + Z.^2);
    
    [ymin, iind] = min(DIST);
    
    mindist(ii) = ymin;
    minind(ii) = iind;
    
end

% sort the solutions

[orb_sort, orb_sort_ind] = sort(mindist);

% plot the top 10?
numplot = 10;

% plot Earth LEO GEO, etc
figure; 
% plot3(1 - mu, 0, 0, '*b')

% --- earth ------
% plot3(1 - mu, 0, 0, '.b')
r_radius = 6371*KILOMETERS/a;
x_offset = 1-mu;
plot_earth(r_radius, x_offset)

% ----earth

hold on
plot3(1 - mu + leo_orb(:,1), leo_orb(:,2), leo_orb(:,3), '--k', 'linewidth', 2)
plot3(1 - mu + geo_orb(:,1), geo_orb(:,2), geo_orb(:,3), '--k', 'linewidth', 2)
plot3(1 - mu + moon_orb(:,1), moon_orb(:,2), moon_orb(:,3), '--b', 'linewidth', 2)
plot3(xnum, ynum, znum,'k-', 'linewidth', 2)
xlabel('$x$ [Au]' , 'Interpreter', 'Latex', 'FontSize', font_size)
ylabel('$y$ [Au]' , 'Interpreter', 'Latex', 'FontSize', font_size)
zlabel('$z$ [Au]' , 'Interpreter', 'Latex', 'FontSize', font_size)
set(gca,'fontsize',font_size-4)
grid on
view([-16, 5])

% orbit 12 hits LEO at 327
for ii = 2:numplot
    
    ind_sol = orb_sort_ind(ii);
    xnumvec = sol{ind_sol}.xnumvec;
   
    plot3(xnumvec(:,1),xnumvec(:,2),xnumvec(:,3), 'linewidth', 0.1)
    % plot the closest approach
    plot3(xnumvec(minind(ind_sol),1),xnumvec(minind(ind_sol),2),...
        xnumvec(minind(ind_sol),3), 'ks')
end
plot3(1 - mu + leo_orb(:,1), leo_orb(:,2), leo_orb(:,3), '--k', 'linewidth', 2)
% plot3(1 - mu, 0, 0, '*b')
plot3(1 - mu + geo_orb(:,1), geo_orb(:,2), geo_orb(:,3), '--k', 'linewidth', 2)
plot3(1 - mu + moon_orb(:,1), moon_orb(:,2), moon_orb(:,3), '--b', 'linewidth', 2)
plot3(xnum, ynum, znum,'k-', 'linewidth', 2)

leg = legend({'Earth', 'LEO', 'GEO','Moon','nonlinear halo'}, 'Interpreter', 'Latex',...
    'location', 'best');
leg.FontSize = font_size;
axis equal
xlim([.992, 1.015]);
ylim([-7e-3 7e-3]);
zlim([-1e-3 1e-3]);
    
%% Test for direct delta v at manifold
close all;
% plot of the results
results = {};

for ii = 1:35
    
    
% non dimensionalized system rotation
omega_rot = [0, 0, 1]';


% non dimensionalized velocity of earth
v_earth_inert = cross(omega_rot, [1-mu, 0, 0]');

% Obtain desired manifold

% manifold solution
mani_ind =ii; %19; %11 + 0*30 + 0*19;
sol_id = orb_sort_ind(mani_ind);

% obtain the state history of that manifold
xvecnum = sol{sol_id}.xnumvec;

% obtain the position closest to earth
x_state_mani = xvecnum(minind(sol_id),1:6)';

% obtain the velocity in rotating frame about sol
v_mani_rot_sol = x_state_mani(4:6);

% obtain rorb (the vector from the center of earth to the manifold)
r_orb = x_state_mani(1:3) - [1-mu; 0; 0];

r_orb_mag = norm(r_orb); % non dimensionalized units

% convert the velocity of the manifold to Earth centered velocity (non dim)
v_mani_inert_sol =  v_mani_rot_sol + cross(omega_rot, x_state_mani(1:3));

% velocity of the manifold w.r.t earth
v_mani_inert_earth = v_mani_inert_sol - v_earth_inert;

% unit vector in the direction of r_orb
e_r = r_orb/norm(r_orb);

% compute the projection vector from the inertial velocity of the manifold
% to the plane with a normal at r
proj_rorb_vmani = v_mani_inert_earth'*e_r * e_r;

% in plane conponent:
inplane_v_mani = v_mani_inert_earth - proj_rorb_vmani;

% unit vector in the direction of the velocity 
e_inplane_v = inplane_v_mani/norm(inplane_v_mani);


% ----- circular orbit calculations -------
% compute the magnitude of circular orbit in m/s at radius r_orb;

v_orb_mag = sqrt(G*M_earth/(r_orb_mag*a)); % m/s 

% non dimensionalization
v_orb_mag = v_orb_mag/a*YEARS/(2*pi); % unitless

v_orb = v_orb_mag*e_inplane_v; % unitless

% ------ converting to rotational baricenter frame ----

v_orb_inert_sol = v_orb + v_earth_inert;

v_orb_rot_sol = v_orb_inert_sol - cross(omega_rot, x_state_mani(1:3));

x_state_orb = [x_state_mani(1:3); v_orb_rot_sol];

opts = odeset('RelTol',1e-14,'AbsTol',1e-14);
% tspan = [0, 4.5*pi];
t_period = 2*pi*sqrt(r_orb_mag^3/mu);

tspan = [0 t_period];

[t,xnumvec_orb] = ode45(@(t,x) eom_3body_rot_fun(t,x,mu), tspan, x_state_orb, opts);


% delta-v

delta_v = v_mani_rot_sol - v_orb_rot_sol;

delta_v_mag = norm(delta_v)*a/YEARS*(2*pi); % m/s

fprintf(['Delta-V to arrive at manifold ' num2str(delta_v_mag/KILOMETERS) ' km/s \n'])  

% -----Total Delta v ----- (hohmann transfer)

% define parking orbit

r1 = (6371+400)*KILOMETERS/a;


% hohman transfer % assumed in the same plane
r2 = r_orb_mag;

% delta_v1 magnitude
delta_v1_mag = sqrt(mu/r1)*(sqrt(2*r2/(r1+r2))-1);


% magnitude of velocity at apogee
v2_orb_mag = sqrt(mu*2*r1/(r2*(r2+r1)));

v2_orb_inert_sol = v2_orb_mag*e_inplane_v + v_earth_inert;

v2_orb_rot_sol = v2_orb_inert_sol - cross(omega_rot, x_state_mani(1:3));

% transfer time
t_transfer = 1/2*2*pi*sqrt((0.5*(r1+r2)*a)^3/(G*M_earth)); % [s]

% initial condition for transfer hohmann orbit
x_state_transorb = [r1*e_r; delta_v1_mag*e_inplane_v];

% delta_v2
delta_v2 = v_mani_rot_sol - v2_orb_rot_sol;

delta_v2_mag = norm(delta_v2);

total_delta_v = (delta_v2_mag + delta_v1_mag)*a/YEARS*2*pi;

fprintf(['Delta-V total ' num2str(total_delta_v/KILOMETERS) ' km/s \n'])  

% --- inclination ----
[z_min, zind_min] = min(xnumvec_orb(:,3));
[z_max, zind_max] = max(xnumvec_orb(:,3));
vec_min = xnumvec_orb(zind_min,1:3)' - [1-mu;0;0];
vec_max = xnumvec_orb(zind_max,1:3)' - [1-mu;0;0];
delvec = vec_max-vec_min;
delvec = delvec/norm(delvec);
inc = acos(delvec'*[0;0;1]);

fprintf(['Orbit inclination ' num2str(90 - inc*180/pi) ' [deg] \n'])  


c3 = 1/2*norm(v_mani_inert_earth*a/YEARS*2*pi)^2 - G*M_earth/(r_orb_mag*a); %m^2/s^2

fprintf(['C3 value ' num2str(c3/KILOMETERS^2) ' [km^2/s^2] \n'])  


% Store all the results
% results nb
results.r_mag_hoh(:,ii) = r2*a/KILOMETERS;
results.mani_id(:,ii) = mani_ind;
results.sol_id(:,ii) = orb_sort_ind(mani_ind);
results.delta_v1_mag(:,ii) = delta_v1_mag;
results.delta_v2_mag(:,ii) = delta_v2_mag;
results.total_delta_v(:,ii) = total_delta_v;
results.c3(:,ii) = c3;
results.inc(:,ii) = pi/2 - inc;
results.x_state_orb(:,ii) = x_state_orb;
results.t_transfer(:,ii) = t_transfer/DAYS;
results.t_manifold(:, ii) = -sol{sol_id}.t(minind(sol_id))/(2*pi)*YEARS/DAYS;
results.t_total(:,ii) = t_transfer/DAYS-sol{sol_id}.t(minind(sol_id))/(2*pi)*YEARS/DAYS;


% plot the solution
if 1
figure(1);
hold on
% --- earth ------
% plot3(1 - mu, 0, 0, '.b')
r_radius = 6371*KILOMETERS/a;
x_offset = (1-mu);
plot_earth(r_radius, x_offset)

% ----earth

plot3(xnumvec_orb(:,1), xnumvec_orb(:,2), xnumvec_orb(:,3), 'k')
%  -----manifold------
ind_sol = orb_sort_ind(mani_ind);
xnumvec_mani = sol{ind_sol}.xnumvec;
plot3(xnumvec_mani(:,1),xnumvec_mani(:,2),xnumvec_mani(:,3), '--b')
% plot the closest approach
plot3(xnumvec_mani(minind(ind_sol),1),xnumvec_mani(minind(ind_sol),2),...
        xnumvec_mani(minind(ind_sol),3), 'k.')
    
% -----manifold end ----
axis equal
xlabel('$x$ [Au]' , 'Interpreter', 'Latex', 'FontSize', font_size)
ylabel('$y$ [Au]' , 'Interpreter', 'Latex', 'FontSize', font_size)
zlabel('$z$ [Au]' , 'Interpreter', 'Latex', 'FontSize', font_size)
set(gca,'fontsize',font_size-4)
grid on

leg = legend({'Earth', 'Parking Orb', 'Mainfold', 'Int. Point'}, 'Interpreter', 'Latex',...
    'location', 'best');
leg.FontSize = font_size;
end

% Earth centered view 5 radia
if 1 
    
num_radia = 5;

figure(2);
hold on
% --- earth ------
% plot3(1 - mu, 0, 0, '.b')
r_radius = 6371*KILOMETERS/a;
x_offset = (1-mu)*0;
plot_earth(1, x_offset)

% ----earth

plot3((xnumvec_orb(:,1)- (1-mu))/r_radius, xnumvec_orb(:,2)/r_radius, xnumvec_orb(:,3)/r_radius, 'k')
%  -----manifold------
ind_sol = orb_sort_ind(mani_ind);
xnumvec_mani = sol{ind_sol}.xnumvec;
plot3((xnumvec_mani(:,1)-(1-mu))/r_radius,xnumvec_mani(:,2)/r_radius,xnumvec_mani(:,3)/r_radius, '--b')
% plot the closest approach
plot3((xnumvec_mani(minind(ind_sol),1)-(1-mu))/r_radius,xnumvec_mani(minind(ind_sol),2)/r_radius,...
        xnumvec_mani(minind(ind_sol),3)/r_radius, 'k.')
    
% -----manifold end ----
axis equal
xlabel('$x$ [$R_E$]' , 'Interpreter', 'Latex', 'FontSize', font_size)
ylabel('$y$ [$R_E$]' , 'Interpreter', 'Latex', 'FontSize', font_size)
zlabel('$z$ [$R_E$]' , 'Interpreter', 'Latex', 'FontSize', font_size)
set(gca,'fontsize',font_size-4)
grid on

xlim([-num_radia, num_radia]);
ylim([-num_radia, num_radia]);
zlim([-num_radia, num_radia]);
leg = legend({'Earth', 'Parking Orb', 'Mainfold', 'Int. Point'}, 'Interpreter', 'Latex');
leg.FontSize = font_size;

end

end

results.sol = sol;
