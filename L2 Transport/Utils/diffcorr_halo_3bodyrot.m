function [x0_dc,t_per, num_iter] = diffcorr_halo_3bodyrot(mu, x0_lin, tol)
%[x0_dc,t_per, num_iter] = DIFFCORR_HALO_3BODYROT(mu, x0_lin, tol)
%Differential corrector for halo orbit in 3 body rotational
%
% Purpose: Obtains the nonlinear initial condition for a halo orbit in a 3
% body rotation problem at the respective tolerance. The tolerance is
% stated as the magnitude of the xy plane crossing velocity.
%
% Inputs:
% mu        [1x1] Gravitational parameter (mass ratio parameter) M2/(M1+M2)
%           given that M1 > M2
% x0_lin    [6x1] Initial solution based on Richardson's 3rd order method
% tol       [1x1] Tolerance for the magnitude of x and z velocity when it
%           crosses the xz plane
%
% Outputs: 
% x0_dc     [6x1] Initial state vector for nonlinear results from
%           differential corrector
% t_per     [1x1] period of the orbit defiend as 2*time of corssing
% num_iter  [1x1] number of iterations needed for the DC to converge to
%           specified tolerance
%
%
% Co-dependencies
%  diffcorr_3body_zx_cross_fun(t, xvec) method that stops ode45 when there
%  is an xz plane crossing
%  eom_3bodystate_STM_fun(t,xvec, mu) main eom function for ode45
%  eom_3body_rot_fun(t,xvec, mu) eom function for 3rd body problem only
%
% Source:
%  [1] Richardson, D.L. Celestial Mechanics (1980) 22: 241. "Analytic
%  Construction of Periodic Orbits about the Col;inear Points."
%  https://doi.org/10.1007/BF01229511
% 
% Author: Alejandro Cabrales H

plot_flag = 0;

if plot_flag
    figure
end

xstate_init = x0_lin;

STM_0 = eye(6);

num_iter = 0;

ErrTol = 1;

while ErrTol > tol && num_iter <= 40

xcombined_init = [xstate_init; STM_0(:)];

opts = odeset('RelTol',1e-13,'AbsTol',1e-14, 'Event', @diffcorr_3body_zx_cross_fun);
tspan = [0, 50];

[t,xcombinednumvec, te,ye, ie] = ode45(@(t,x) eom_3bodystate_STM_fun(t,x,mu), tspan, xcombined_init, opts);

% Unpack end states
x0 = xcombinednumvec(1, 1:6)';
xf = xcombinednumvec(end, 1:6)';

% Compute dXvec/xt at xf
xdot = eom_3body_rot_fun(0, xf, mu);

phi_tf_t0 = reshape(xcombinednumvec(end, 7:end), 6,6);

% comptue deltax = phi*x0 + dx(tf)/xt * d(T/2);

% set dz = 0
phi_3d = [phi_tf_t0(4,1) - xdot(4)/xf(5)*phi_tf_t0(2,1), phi_tf_t0(4,5) - xdot(4)/xf(5)*phi_tf_t0(2,5);
    phi_tf_t0(6,1) - xdot(6)/xf(5)*phi_tf_t0(2,1), phi_tf_t0(6,5) - xdot(6)/xf(5)*phi_tf_t0(2,5)];


if rank(phi_3d) == 1
    % utilize different correction in which dx = 0
    
    % set dx = 0
    phi_3d = [phi_tf_t0(4,3) - xdot(4)/xf(5)*phi_tf_t0(2,3), phi_tf_t0(4,5) - xdot(4)/xf(5)*phi_tf_t0(2,5);
        phi_tf_t0(6,3) - xdot(6)/xf(5)*phi_tf_t0(2,3), phi_tf_t0(6,5) - xdot(6)/xf(5)*phi_tf_t0(2,5)];
    
    % Delta Vector
    delta = (phi_3d)\[-xf(4); -xf(6)];
    
    xstate_init(3) = xstate_init(3) + delta(1);
    xstate_init(5) = xstate_init(5) + delta(2);
    
else % matrix is full rank so delta will have a solution
    
    delta = (phi_3d)\[-xf(4); -xf(6)];
    
    xstate_init(1) = xstate_init(1) + delta(1);
    xstate_init(5) = xstate_init(5) + delta(2);
end


if plot_flag
font_size = 14;
xnum = xcombinednumvec(:,1);
ynum = xcombinednumvec(:,2);
znum = xcombinednumvec(:,3);
    
% Y Z Axis
subplot(2,2,1)
plot(ynum, znum,'k-')
hold on
% plot(y,z, 'b--')
xlabel('$y$ [$\gamma_L$]' , 'Interpreter', 'Latex', 'FontSize', font_size)
ylabel('$z$ [$\gamma_L$]' , 'Interpreter', 'Latex', 'FontSize', font_size)
set(gca,'fontsize',font_size-4)
grid on
% axis equal
% X Z axis
subplot(2,2,2)
plot(xnum, znum,'k-')
hold on
% plot(x, z,'b--')
xlabel('$x$ [$\gamma_L$]' , 'Interpreter', 'Latex', 'FontSize', font_size)
ylabel('$z$ [$\gamma_L$]' , 'Interpreter', 'Latex', 'FontSize', font_size)
set(gca,'fontsize',font_size-4)
grid on
% axis equal

% X Y axis
subplot(2,2,3)
plot(xnum, ynum,'k-')
hold on
% plot(x, y,'b--')
xlabel('$x$ [$\gamma_L$]' , 'Interpreter', 'Latex', 'FontSize', font_size)
ylabel('$y$ [$\gamma_L$]' , 'Interpreter', 'Latex', 'FontSize', font_size)
set(gca,'fontsize',font_size-4)
grid on
% axis equal
xlim([-1 1]);
ylim([-1 1]);

% 3d
subplot(2,2,4)
plot3(xnum, ynum, znum,'k-')
hold on
% plot3(x, y, z,'b--')
xlabel('$x$ [$\gamma_L$]' , 'Interpreter', 'Latex', 'FontSize', font_size)
ylabel('$y$ [$\gamma_L$]' , 'Interpreter', 'Latex', 'FontSize', font_size)
zlabel('$z$ [$\gamma_L$]' , 'Interpreter', 'Latex', 'FontSize', font_size)
set(gca,'fontsize',font_size-4)
grid on
% axis equal
view([-41 21])
drawnow
end


num_iter = num_iter +1;
ErrTol = sqrt(abs(xf(4))^2 + abs(xf(6))^2);
end

x0_dc = xstate_init;

t_per = 2*max(te);

end

