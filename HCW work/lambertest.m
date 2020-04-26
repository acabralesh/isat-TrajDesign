% lambertest

close all;
% addpath('Utils')
% addpath('../Units')
font_size = 22;
% universal values

mu = 3.986004418*10^(14)*METERS^3/SECONDS^2;

earth_radius = 6371*KILOMETERS;

v_earth = 0*491*METERS/SECONDS;

opts = odeset('RelTol',1e-10,'AbsTol',1e-10);

%% Satellite orbits

% Obtain orbit for sat 1
r_orb1 = earth_radius + 1000*KILOMETERS;

v_orb1 = sqrt(mu/r_orb1);

e_r1 = [1 0 0]';

d1 = [0 1 0]';

xinit1 = [r_orb1*e_r1; v_orb1*d1];

T1 = 2*pi*sqrt(r_orb1^3/mu);

t1 = [0 T1];

% Obtain orbit for sat 2
r_orb2 = earth_radius + 2000*KILOMETERS;

v_orb2 = sqrt(mu/r_orb2);

e_r2 = [0 1 0]'; e_r2 = e_r2/norm(e_r2);

d2 = [-1 0 0]';

xinit2 = [r_orb2*e_r2; v_orb2*d2];

T2 = 2*pi*sqrt(r_orb2^3/mu);

t2 = [0 T2];


[t1,xvec_sc1] = ode45(@(t,x) odekep(t,x,0,mu,1), ...
        t1, xinit1, opts);
    
[t2,xvec_sc2] = ode45(@(t,x) odekep(t,x,0,mu,1), ...
        t2, xinit2, opts);
    
%% Manifold test
  
  load('xmanifold.mat')    
t = xmanifold.t;
xinit_mani = xmanifold.xnumvec_eci_m(1:6,:);
%% Lamberts Test

t_ind = find(t2 > T2/2,1);
r1_vec = xvec_sc1(1,1:3)';
v1_vec = xvec_sc1(1,4:6)';
x1_state = [r1_vec; v1_vec];

r2_vec = xvec_sc2(t_ind,1:3)';
v2_vec = xvec_sc2(t_ind,4:6)';
x2_state = [r2_vec; v2_vec];

% minimum a lambert problem solution 
[v1_lamvec,v2_lamvec, T_lamb] = lambert_mina(mu,x1_state, x2_state);

% general lambert problem solution
% [vi, vf] = glambert(mu, x1_state, x2_state, 5*T_lamb, -1);   

xinit_lamb = [r1_vec; v1_lamvec];
% xinit_lamb = [r1_vec; vi];

t_lamb = [0 T_lamb];

%% Integrate
    
[t_lamb,xvec_lam] = ode45(@(t,x) odekep(t,x,0,mu,1), ...
        t_lamb, xinit_lamb, opts);

%% Plotting    
figure;

plot_earth(1, 0);
hold on
plot3(xvec_sc1(:,1)/earth_radius, xvec_sc1(:,2)/earth_radius, xvec_sc1(:,3)/earth_radius,'k--', 'linewidth', 1)
plot3(xvec_sc2(:,1)/earth_radius, xvec_sc2(:,2)/earth_radius, xvec_sc2(:,3)/earth_radius,'k-', 'linewidth', 1)

% lamb
plot3(r1_vec(1)/earth_radius, r1_vec(2)/earth_radius, r1_vec(3)/earth_radius,'^k', 'linewidth', 1)
plot3(r2_vec(1)/earth_radius, r2_vec(2)/earth_radius, r2_vec(3)/earth_radius,'sk', 'linewidth', 1)

plot3(xvec_lam(:,1)/earth_radius, xvec_lam(:,2)/earth_radius, xvec_lam(:,3)/earth_radius,'k-', 'linewidth', 1)

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