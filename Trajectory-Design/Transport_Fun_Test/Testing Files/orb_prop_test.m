% Testing propagator for CR3BP
% 
load('EML1_manifolds.mat')
idx = 1


xvecnum_test = EML1_manifolds.sol{idx}.xnumvec;
a = EML1_manifolds.a;
x_state_orb = xvecnum_test(end,:);
mu = EML1_manifolds.mu;
tend = -EML1_manifolds.sol{idx}.tnumvec(end);


pos_units_convert = a;
vel_units_convert = a/(27.3*DAYS)*2*pi;
time_units_convert = 2*pi/(27.3*DAYS);


opts = odeset('RelTol',1e-14,'AbsTol',1e-14);
% tspan = [0, 4.5*pi];
t_period =tend;

tspan = [0 t_period];

[t,xnumvec_orb] = ode45(@(t,x) eom_3body_rot_fun(t,x,mu), tspan, x_state_orb, opts);

figure;
plot3(xnumvec_orb(:,1)*a/earthRadius, xnumvec_orb(:,2)*a/earthRadius, xnumvec_orb(:,3)*a/earthRadius, 'k-', 'linewidt', 1)
hold on
plot3(xnumvec_orb(1,1)*a/earthRadius, xnumvec_orb(1,2)*a/earthRadius, xnumvec_orb(1,3)*a/earthRadius, 'k*')
plot_earth(1, -mu*a/earthRadius)

%% Test that lamberts work

% target position
target_pos = [-earthRadius-mu*a; 0; 0];
target_vel = [0 -410.9632 0]';

% ---------------


lei_R_cei = eye(3);

taget_loc = xnumvec_orb(1,:)';
omega_n = [0 0 1]';
[sc_pos_lei,sc_vel_lei] = ber_to_bei(eye(3), omega_n,taget_loc(1:3), taget_loc(4:6));

% Convert LER to CEI - - Earth is the initial center
r_earth_lei = eye(3)*[-mu;0;0];
        
[sc_pos_cei,sc_vel_cei] = bei_to_cei(lei_R_cei, omega_n,r_earth_lei, sc_pos_lei, sc_vel_lei);

x2_state = [sc_pos_cei*pos_units_convert; sc_vel_cei*vel_units_convert];

x1_state = [target_pos; target_vel];
% minimum a lambert problem solution
Gmu = 3.986004418*10^(14)*METERS^3/SECONDS^2;
[v1_lamvec,v2_lamvec, T_lamb] = lambert_mina(Gmu, x1_state, x2_state);

% Convert to CER again
 % non dimensionalize output
 target_pos_nd = target_pos/pos_units_convert;
 v1_lamvec_nd = (v1_lamvec)/vel_units_convert;
 T_lamb_nd = T_lamb*time_units_convert;


[sc_pos_bei,sc_vel_bei] = cei_to_bei(lei_R_cei, omega_n,r_earth_lei, target_pos_nd(:), v1_lamvec_nd(:));

[sc_pos_ber,sc_vel_ber] = bei_to_ber(lei_R_cei, omega_n, sc_pos_bei, sc_vel_bei);


x_state_trans = [sc_pos_ber;sc_vel_ber];
tspan = [0 T_lamb_nd];
[t,xnumvec_trans] = ode45(@(t,x) eom_3body_rot_fun(t,x,mu), tspan, x_state_trans, opts);

[t_trans, xvec_trans] = propagate_lambert_fun(Gmu, T_lamb, [target_pos; v1_lamvec*0.997], opts);

% plot3(xvec_trans(1,1)/earthRadius, xvec_trans(1,2)/earthRadius, xvec_trans(1,3)/earthRadius, 'r*', 'linewidt', 1)
plot3(xvec_trans(:,1)/earthRadius, xvec_trans(:,2)/earthRadius, xvec_trans(:,3)/earthRadius, 'k-', 'linewidt', 1)

plot3(taget_loc(1)*a/earthRadius, taget_loc(2)*a/earthRadius, taget_loc(3)*a/earthRadius, 'ok')
        
 xlabel('$x \ [R_E]$','Interpreter','latex', 'FontSize',20)
 
 ylabel('$y \ [R_E]$','Interpreter','latex', 'FontSize',20)
 
 zlabel('$z \ [R_E]$','Interpreter','latex', 'FontSize',20)
 axis equal
%  xlim([-3 3])
%  ylim([-3 3])
%  zlim([-3 3])
function [t_trans, xvec_trans] = propagate_lambert_fun(mu, T_lamb, xinit_lamb, opts)

[t_trans,xvec_trans] = ode45(@(t,x) odekep(t,x,0,mu,1), ...
    [0 T_lamb], xinit_lamb, opts);

end
