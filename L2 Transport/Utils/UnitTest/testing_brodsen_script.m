% Testing Brosden Script
addpath('Units')
addpath('Utils')
close all
font_size = 18;

G = 6.67408e-11 * METERS^3 * KILOGRAMS^-1 * SECONDS^-2;
M_sun = 1.9884754e30 * KILOGRAMS;
M_earth = 5.972358e24 * KILOGRAMS;
M_moon = 7.34767309e22 * KILOGRAMS;

% Sun - Earth
a = 1.4959787e8 * KILOMETERS;
% Earth - Moon
% a = 384400 * KILOMETERS;

% Sun
mu = M_earth/(M_sun + M_earth);
% mu = 3.04036e-6;

% Earth
% mu = M_moon/(M_earth + M_moon);


% Numerically find lagrange points
options = optimset('TolX', 1e-16);

% Sun
[x_l1] = fzero(@(x) diff_potential_3body_fun(mu, [x,0,0]', 1), .999, options)
% Moon
% [x_l1] = fzero(@(x) diff_potential_3body_fun(mu, [x,0,0]', 1), .7, options);


[x_l2] = fzero(@(x) diff_potential_3body_fun(mu, [x,0,0]', 1), 1.2, options);

% check that the potentials are close to zero
dU_l1 = diff_potential_3body_fun(mu, [x_l1,0,0]', 1:3);

dU_l2 = diff_potential_3body_fun(mu, [x_l2,0,0]', 1:3);

fprintf(['Mass ratio (mu) ' num2str(mu) '\n'])
fprintf(['Location of L1 (non dimen) ' num2str(x_l1) '\n'])
fprintf(['Location of L2 (non dimen) ' num2str(x_l2) '\n'])
fprintf(['Location of L1 from Secondary (Km) ' num2str((x_l1-1)*a/KILOMETERS) '\n'])
fprintf(['Location of L2 from Secondary (Km) ' num2str((x_l2-1)*a/KILOMETERS) '\n'])

% modified
L1 = roots([1,-(3-mu),(3-2*mu),-mu,2*mu,-mu]);
L2 = roots([1,(3-mu),(3-2*mu),-mu,-2*mu,-mu]);
L3 = roots([1,(2+mu),(1+2*mu),-(1-mu),-2*(1-mu),-(1-mu)]);

%%% x-position of Lagrange points in rotating frame
xL1 = L1(find(imag(L1)==0 & real(L1)>0));
xL2 = L2(find(imag(L2)==0 & real(L2)>0));
xL3 = L3(find(imag(L3)==0 & real(L3)>0));

%%% distance from Sun (mass 1 - primary) to Lagrange points 1, 2, and 3
r1_L1 = abs(1-xL1);
r1_L2 = abs(1+xL2);
r1_L3 = abs(xL3);


x_l1 = r1_L1
%modified

%% Contour Plot 
% Sun
xvec = linspace(x_l1-.01, x_l2+.01, 100);
yvec = linspace(-.01, .01, 100);

% % Moon
% xvec = linspace(x_l1-.1, x_l2+.1, 100);
% yvec = linspace(-.1, .1, 100);

[X, Y] = meshgrid(xvec, yvec); 

Z = zeros(size(X));

U = potential_3body_fun(mu, X, Y, Z);

figure; 
contour(X,Y, U, linspace(1.5003, 1.5007, 20)); % Sun - Earth 
% contour(X,Y, U, linspace(1.500, 1.7, 20)); % Earth - Moon 
hold on
plot(x_l1, 0, '*k')
plot(x_l2, 0, '*k')
xlabel('$x$ (non-dimensional)' , 'Interpreter', 'Latex', 'FontSize', font_size)
ylabel('$y$ (non-dimensional)' , 'Interpreter', 'Latex', 'FontSize', font_size)
title({'Zero Velocity Contour Plots'; ['$\mu =$', num2str(mu)]} , 'Interpreter', 'Latex', 'FontSize', font_size)
set(gca,'fontsize',font_size-4)
grid on

%% First test of Lagrange point orbit
% initial state

if 0 
xinit = [0.988837156239640; 0; 0.000810869832364; 0; 0.008939405772417; 0];

opts = odeset('RelTol',1e-12,'AbsTol',1e-12);
tspan = [0, 3];

[t,x] = ode45(@(t,x) eom_3body_rot_fun(t,x,mu), tspan, xinit, opts);

figure;
plot3(x(:,1), x(:,2), x(:,3)); hold on; plot3(x(1,1), x(1,2), x(1,3), '*r')
grid on
xlabel('x')
ylabel('y')
zlabel('z')

end
%% Richardson Testing

% -----inputs-----
% deffined values
% mu = 3.040357143e-6; % gravitational parametesr
% gamma_l = 1.001090475e-2; % ratio of distance from l1 to l2
% a = 1.495978714e8*KILOMETERS;
% 
% mu = 3.04036e-6;
% gamma_l = 1.00109e-2;
% a = 1.49598e8*KILOMETERS;



gamma_l = (1-x_l1);
% x_l1 = 1-gamma_l;
% gamma_l = (x_l2-1);
libind = 1; % wether its libration point 1 or 2
Az = 110000*KILOMETERS/(gamma_l*a); % amplitude Az desired

% unsure about phi
phi = 0;

n = 1; % 1 for class 1, 3 for class 2

% ----- Begin Function -----

% c coefficients
c2 = c_libration(mu, gamma_l, libind, 2);
c3 = c_libration(mu, gamma_l, libind, 3);
c4  = c_libration(mu, gamma_l, libind, 4);

% eigen values
% obtain the solution to lamb^4 + (2-c_2)*lamb^2 + (1+2*c_2)*(1-c_2) = 0
eig_l = roots([1 0 (2-c2) 0 (1+2*c2)*(1-c2)]);

lamb = max(imag(eig_l));

% correction coefficient for y
k = 1/(2*lamb)*(lamb^2 + 1 + 2*c2);

% Correction for conmesurate condition on z
delta = lamb^2 - c2;

% Coefficients for Richardson 3rd order

% d1 d2 coefficients
d1 = 3*lamb^2/k*(k*(6*lamb^2 - 1) - 2*lamb);

d2 = 8*lamb^2/k*(k*(11*lamb^2 - 1) - 2*lamb);

% a21 a22 a23 a24 coefficients
a21 = 3*c3*(k^2 - 2)/(4*(1 + 2*c2));

a22 = 3*c3/(4*(1 + 2*c2));

a23 = -(3*c3*lamb/(4*k*d1))*(3*k^3*lamb - 6*k*(k - lamb) + 4);

a24 = -(3*c3*lamb/(4*k*d1))*(2 + 3*k*lamb);

% d21 d31 d32 coefficients 
d21 = -c3/(2*lamb^2);

d31 = 3/(64*lamb^2)*(4*c3*a24 + c4);

d32 = 3/(64*lamb^2)*(4*c3*(a23 - d21) + c4*(4 + k^2));

% b coefficients pt 1
b21 = -(3*c3*lamb/(2*d1))*(3*k*lamb - 4);

b22 = 3*c3*lamb/d1;

b31 = (3/(8*d2))*( 8*lamb*( 3*c3*(k*b21 - 2*a23) - c4*(2 + 3*k^2)) + ...
    (9*lamb^2 + 1 + 2*c2)*( 4*c3*(k*a23 - b21) + k*c4*(4 + k^2)) );

b32 = 1/d2*( 9*lamb*( c3*(k*b22 + d21 - 2*a24) - c4) + ...
    3/8*(9*lamb^2 + 1 + 2*c2)*( 4*c3*(k*a24 - b22) + k*c4) );

% a1 a2 coefficients
a1 = -3/2*c3*(2*a21 + a23 + 5*d21) - 3/8*c4*(12 - k^2);

a2 = 3/2*c3*(a24 - 2*a22) + 9/8*c4;

% a31 a32 coefficients
a31 = -(9*lamb/(4*d2))*( 4*c3*(k*a23 - b21) + k*c4*(4 + k^2) ) + ...
    ((9*lamb^2 + 1 - c2)/(2*d2))*( 3*c3*(2*a23 - k*b21) + c4*(2 + 3*k^2));

a32 = -1/d2*( 9*lamb/4*( 4*c3*(k*a24 - b22) + k*c4 ) + ...
    3/2*(9*lamb^2 + 1 - c2)*( c3*(k*b22 + d21 - 2*a24) - c4 ) );

% s1 s2 coefficients
s1 = (1/(2*lamb*( lamb*(1 + k^2) - 2*k)))*( 3/2*c3*( 2*a21*(k^2 - 2) - ...
    a23*(k^2 + 2) - 2*k*b21) - 3/8*c4*(3*k^4 - 8*k^2 + 8));

s2 = (1/(2*lamb*( lamb*(1 + k^2) - 2*k)))*( 3/2*c3*( 2*a22*(k^2 - 2) + ...
    a24*(k^2 + 2) + 2*k*b22 + 5*d21 ) + 3/8*c4*(12 - k^2) );

% l1 l2 coefficients
l1 = a1 + 2*lamb^2*s1;

l2 = a2 + 2*lamb^2*s2;

% Generate the amplitude 
Ax = sqrt((-l2*Az^2 - delta)/l1);

% w1 w2
w1 = 0;

w2 = s1*Ax^2 + s2*Az^2;

w = 1+ w1 + w2; 

% phi, psi
psi = phi + n*pi/2;

% Plotting function

tau = w*linspace(0, 3.1, 100);

tau1 = lamb*tau + phi;

deltan = 2 - n;

x = a21*Ax^2 + a22*Az^2 - Ax*cos(tau1) + (a23*Ax^2 - a24*Az^2)*cos(2*tau1) + ...
    (a31*Ax^3 - a32*Ax*Az^2)*cos(3*tau1);

y = k*Ax*sin(tau1) + (b21*Ax^2 - b22*Az^2)*sin(2*tau1) + ...
    (b31*Ax^3 - b32*Ax*Az^2)*sin(3*tau1);

z = deltan*Az*cos(tau1) + deltan*d21*Ax*Az*(cos(2*tau1) - 3) + ...
    deltan*(d32*Az*Ax^2 - d31*Az^3)*cos(3*tau1);

xd = w*lamb*Ax*sin(tau1) - 2*w*lamb*(a23*Ax^2 - a24*Az^2)*sin(2*tau1) - ...
    3*w*lamb*(a31*Ax^3 - a32*Ax*Az^2)*sin(3*tau1);

yd = w*lamb*k*Ax*cos(tau1) + 2*w*lamb*(b21*Ax^2 - b22*Az^2)*cos(2*tau1) + ...
    3*lamb*(b31*Ax^3 - b32*Ax*Az^2)*cos(3*tau1);

zd = -w*lamb*deltan*Az*sin(tau1) - 2*w*lamb*deltan*d21*Ax*Az*sin(2*tau1) - ... 
    3*w*lamb*deltan*(d32*Az*Ax^2 - d31*Az^3)*sin(3*tau1);

% plotting

% converting to KM
% x = x*gamma_l*a/KILOMETERS;
% y = y*gamma_l*a/KILOMETERS;
% z = z*gamma_l*a/KILOMETERS;

figure
% Y Z Axis
subplot(2,2,1)
plot(y, z,'k-')
xlabel('$y$ [Km]' , 'Interpreter', 'Latex', 'FontSize', font_size)
ylabel('$z$ [Km]' , 'Interpreter', 'Latex', 'FontSize', font_size)
set(gca,'fontsize',font_size-4)
grid on
axis equal
% X Z axis
subplot(2,2,2)
plot(x, z,'k-')
xlabel('$x$ [Km]' , 'Interpreter', 'Latex', 'FontSize', font_size)
ylabel('$z$ [Km]' , 'Interpreter', 'Latex', 'FontSize', font_size)
set(gca,'fontsize',font_size-4)
grid on
axis equal

% X Y axis
subplot(2,2,3)
plot(x, y,'k-')
xlabel('$x$ [Km]' , 'Interpreter', 'Latex', 'FontSize', font_size)
ylabel('$y$ [Km]' , 'Interpreter', 'Latex', 'FontSize', font_size)
set(gca,'fontsize',font_size-4)
grid on
axis equal
% 3d
subplot(2,2,4)
plot3(x, y, z,'k-')
xlabel('$x$ [Km]' , 'Interpreter', 'Latex', 'FontSize', font_size)
ylabel('$y$ [Km]' , 'Interpreter', 'Latex', 'FontSize', font_size)
zlabel('$z$ [Km]' , 'Interpreter', 'Latex', 'FontSize', font_size)
set(gca,'fontsize',font_size-4)
grid on
axis equal
view([-41 21])

if 0
    
    % converting to KM
x = x*gamma_l*a/KILOMETERS;
y = y*gamma_l*a/KILOMETERS;
z = z*gamma_l*a/KILOMETERS;
ylimbnds = [min(y) max(y)] *1.1;
xlimbnds = [min(x) max(x)] *1.1;
zlimbnds = [min(z) max(z)] *1.1;
   
figure; plot3(x, y, z,'k-', 'linewidth', 1); hold on; 
plot3(x, y*0 + ones(size(z))*ylimbnds(2), z,'--', 'Color', [0   0.447000000000000   0.741000000000000]);
plot3(x, y, z*0 + ones(size(z))*zlimbnds(1),'--', 'Color', [0   0.447000000000000   0.741000000000000]);
plot3(x*0 + ones(size(z))*xlimbnds(2), y, z,'--', 'Color', [0   0.447000000000000   0.741000000000000]);
 zlim(zlimbnds);xlim(xlimbnds); ylim(ylimbnds)
xlabel('$x$ [Km]' , 'Interpreter', 'Latex', 'FontSize', font_size)
ylabel('$y$ [Km]' , 'Interpreter', 'Latex', 'FontSize', font_size)
zlabel('$z$ [Km]' , 'Interpreter', 'Latex', 'FontSize', font_size)
set(gca,'fontsize',font_size-4)
grid on
axis equal
    
end

%% testing 3rd order approx
(x_l1-1)*a
% have xvec = xvec(tau1 = 0)
xinit = [x(1); y(1); z(1); xd(1); yd(1); zd(1)]*gamma_l + [1-mu-gamma_l 0 0 0 0 0]';

xinit = ([x(1); y(1); z(1); xd(1); yd(1); zd(1)]*(1-x_l1)*a + [1-mu-gamma_l 0 0 0 0 0]'*a)/a;


opts = odeset('RelTol',1e-12,'AbsTol',1e-12);
tspan = [0, 3];

[t,xnumvec] = ode45(@(t,x) eom_3body_rot_fun(t,x,mu), tspan, xinit, opts);

xnum = (xnumvec(:,1)-x_l1)/gamma_l;
ynum = xnumvec(:,2)/gamma_l;
znum = xnumvec(:,3)/gamma_l;

figure
% Y Z Axis
subplot(2,2,1)
plot(ynum, znum,'k-')
hold on
plot(y,z, 'b--')
xlabel('$y$ [Km]' , 'Interpreter', 'Latex', 'FontSize', font_size)
ylabel('$z$ [Km]' , 'Interpreter', 'Latex', 'FontSize', font_size)
set(gca,'fontsize',font_size-4)
grid on
axis equal
% X Z axis
subplot(2,2,2)
plot(xnum, znum,'k-')
hold on
plot(x, z,'b--')
xlabel('$x$ [Km]' , 'Interpreter', 'Latex', 'FontSize', font_size)
ylabel('$z$ [Km]' , 'Interpreter', 'Latex', 'FontSize', font_size)
set(gca,'fontsize',font_size-4)
grid on
axis equal

% X Y axis
subplot(2,2,3)
plot(xnum, ynum,'k-')
hold on
plot(x, y,'b--')
xlabel('$x$ [Km]' , 'Interpreter', 'Latex', 'FontSize', font_size)
ylabel('$y$ [Km]' , 'Interpreter', 'Latex', 'FontSize', font_size)
set(gca,'fontsize',font_size-4)
grid on
axis equal
% 3d
subplot(2,2,4)
plot3(xnum, ynum, znum,'k-')
hold on
plot3(x, y, z,'b--')
xlabel('$x$ [Km]' , 'Interpreter', 'Latex', 'FontSize', font_size)
ylabel('$y$ [Km]' , 'Interpreter', 'Latex', 'FontSize', font_size)
zlabel('$z$ [Km]' , 'Interpreter', 'Latex', 'FontSize', font_size)
set(gca,'fontsize',font_size-4)
grid on
axis equal
view([-41 21])

%% Testing State Transition for differential corrector


xstate_init = xinit;

STM_0 = eye(6);
figure
for ii = 1:1e4

xcombined_init = [xstate_init; STM_0(:)];

opts = odeset('RelTol',1e-14,'AbsTol',1e-15, 'Event', @diffcorr_3body_zx_cross_fun);
tspan = [0, 50*1.50623074975922];

[t,xcombinednumvec, te,ye, ie] = ode45(@(t,x) eom_3bodystate_STM_fun(t,x,mu), tspan, xcombined_init, opts);

xnum = (xcombinednumvec(:,1) - x_l1)/gamma_l;
ynum = xcombinednumvec(:,2)/gamma_l;
znum = xcombinednumvec(:,3)/gamma_l;

% xnum = (xcombinednumvec(:,1) );
% ynum = xcombinednumvec(:,2);
% znum = xcombinednumvec(:,3);

% obtain phi
x0 = xcombinednumvec(1, 1:6)'
xf = xcombinednumvec(end, 1:6)'
xdot = eom_3body_rot_fun(0, xf, mu);

phi_tf_t0 = reshape(xcombinednumvec(end, 7:end), 6,6);
% 
% phi_3d = [phi_tf_t0(4,1) - xdot(4)/xf(5)*phi_tf_t0(2,1), phi_tf_t0(4,3) - ...
%     xdot(4)/xf(5)*phi_tf_t0(2,3), phi_tf_t0(4,5) - xdot(4)/xf(5)*phi_tf_t0(2,5);
%     phi_tf_t0(6,1) - xdot(6)/xf(5)*phi_tf_t0(2,1), phi_tf_t0(6,3) - ...
%     xdot(6)/xf(5)*phi_tf_t0(2,3), phi_tf_t0(6,5) - xdot(6)/xf(5)*phi_tf_t0(2,5)];
% 
% delta = pinv(phi_3d)*[-xf(4); -xf(6)];
% 
% 
% % % set dz = 0
% phi_3d = [phi_tf_t0(4,1) - xdot(4)/xf(5)*phi_tf_t0(2,1), phi_tf_t0(4,5) - xdot(4)/xf(5)*phi_tf_t0(2,5);
%     phi_tf_t0(6,1) - xdot(6)/xf(5)*phi_tf_t0(2,1), phi_tf_t0(6,5) - xdot(6)/xf(5)*phi_tf_t0(2,5)];
% 
% delta = inv(phi_3d)*[-xf(4); -xf(6)]

% % set dx = 0
% phi_3d = [phi_tf_t0(4,3) - xdot(4)/xf(5)*phi_tf_t0(2,3), phi_tf_t0(4,5) - xdot(4)/xf(5)*phi_tf_t0(2,5);
%     phi_tf_t0(6,3) - xdot(6)/xf(5)*phi_tf_t0(2,3), phi_tf_t0(6,5) - xdot(6)/xf(5)*phi_tf_t0(2,5)];
% 
% delta = inv(phi_3d)*[-xf(4); -xf(6)]

% % 
phi_wills = [phi_tf_t0(2,1), phi_tf_t0(2,5), xf(5); 
    phi_tf_t0(4,1), phi_tf_t0(4,5), xdot(4);
    phi_tf_t0(6,1), phi_tf_t0(6,5), xdot(6)];

deltawill = inv(phi_wills)*[0; xdot(4); xdot(6)]
% 
% % % 2D solution
% deltay = -xf(4)/( phi_tf_t0(4,5) - xdot(4)/xf(5)*phi_tf_t0(2,5))
% 

% non normalized values
% xnum = (xcombinednumvec(:,1) );
% ynum = xcombinednumvec(:,2);
% znum = xcombinednumvec(:,3);

% opts = odeset('RelTol',1e-14,'AbsTol',1e-15);
% tspan = [0, 3];
% % 
% [t,xcombinednumvec] = ode45(@(t,x) eom_3bodystate_STM_fun(t,x,mu), tspan, xcombined_init, opts);

xnum = (xcombinednumvec(:,1) - x_l1)/gamma_l;
ynum = xcombinednumvec(:,2)/gamma_l;
znum = xcombinednumvec(:,3)/gamma_l;
% xstate_init = xinit;



% xstate_init(1) = xstate_init(1) + delta(1);
% xstate_init(3) = xstate_init(3) + delta(2);
% xstate_init(5) = xstate_init(5) + delta(3);

% xstate_init(1) = xstate_init(1) + delta(1);
% xstate_init(3) = xstate_init(3) + delta(2);
% xstate_init(5) = xstate_init(5) + delta(2);

xstate_init(1) = xstate_init(1) + deltawill(1);
xstate_init(5) = xstate_init(5) + deltawill(2);
% xstate_init(5) = xstate_init(5) + deltay;

if 1 
% Y Z Axis
subplot(2,2,1)
plot(ynum, znum,'k-')
hold on
plot(y,z, 'b--')
xlabel('$y$' , 'Interpreter', 'Latex', 'FontSize', font_size)
ylabel('$z$' , 'Interpreter', 'Latex', 'FontSize', font_size)
set(gca,'fontsize',font_size-4)
grid on
axis equal
% X Z axis
subplot(2,2,2)
plot(xnum, znum,'k-')
hold on
plot(x, z,'b--')
xlabel('$x$' , 'Interpreter', 'Latex', 'FontSize', font_size)
ylabel('$z$' , 'Interpreter', 'Latex', 'FontSize', font_size)
set(gca,'fontsize',font_size-4)
grid on
axis equal

% X Y axis
subplot(2,2,3)
plot(xnum, ynum,'k-')
hold on
plot(x, y,'b--')
xlabel('$x$' , 'Interpreter', 'Latex', 'FontSize', font_size)
ylabel('$y$' , 'Interpreter', 'Latex', 'FontSize', font_size)
set(gca,'fontsize',font_size-4)
grid on
axis equal
xlim([-1 1]);
ylim([-1 1]);

% 3d
subplot(2,2,4)
plot3(xnum, ynum, znum,'k-')
hold on
plot3(x, y, z,'b--')
xlabel('$x$' , 'Interpreter', 'Latex', 'FontSize', font_size)
ylabel('$y$' , 'Interpreter', 'Latex', 'FontSize', font_size)
zlabel('$z$' , 'Interpreter', 'Latex', 'FontSize', font_size)
set(gca,'fontsize',font_size-4)
grid on
axis equal
view([-41 21])
end
drawnow

if abs(xf(4)) < 1e-8 && abs(xf(6)) < 1e-8
    display('completed')
   break 
end



end

%%
opts = odeset('RelTol',1e-14,'AbsTol',1e-15);
tspan = [0, 7];

[t,xcombinednumvec] = ode45(@(t,x) eom_3bodystate_STM_fun(t,x,mu), tspan, xcombined_init, opts);
xnum = (xcombinednumvec(:,1) - x_l1)/gamma_l;
ynum = xcombinednumvec(:,2)/gamma_l;
znum = xcombinednumvec(:,3)/gamma_l;

figure
% Y Z Axis
subplot(2,2,1)
plot(ynum, znum,'k-')
hold on
plot(y,z, 'b--')
xlabel('$y$ [Km]' , 'Interpreter', 'Latex', 'FontSize', font_size)
ylabel('$z$ [Km]' , 'Interpreter', 'Latex', 'FontSize', font_size)
set(gca,'fontsize',font_size-4)
grid on
axis equal
% X Z axis
subplot(2,2,2)
plot(xnum, znum,'k-')
hold on
plot(x, z,'b--')
xlabel('$x$ [Km]' , 'Interpreter', 'Latex', 'FontSize', font_size)
ylabel('$z$ [Km]' , 'Interpreter', 'Latex', 'FontSize', font_size)
set(gca,'fontsize',font_size-4)
grid on
axis equal

% X Y axis
subplot(2,2,3)
plot(xnum, ynum,'k-')
hold on
plot(x, y,'b--')
xlabel('$x$ [Km]' , 'Interpreter', 'Latex', 'FontSize', font_size)
ylabel('$y$ [Km]' , 'Interpreter', 'Latex', 'FontSize', font_size)
set(gca,'fontsize',font_size-4)
grid on
axis equal
% 3d
subplot(2,2,4)
plot3(xnum, ynum, znum,'k-')
hold on
plot3(x, y, z,'b--')
xlabel('$x$ [Km]' , 'Interpreter', 'Latex', 'FontSize', font_size)
ylabel('$y$ [Km]' , 'Interpreter', 'Latex', 'FontSize', font_size)
zlabel('$z$ [Km]' , 'Interpreter', 'Latex', 'FontSize', font_size)
set(gca,'fontsize',font_size-4)
grid on
axis equal
view([-41 21])