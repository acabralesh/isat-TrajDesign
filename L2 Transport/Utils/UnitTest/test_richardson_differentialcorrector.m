%% Testing Differential Corrector
% clc; 
% close all; 
addpath('Utils')
addpath('Units')
font_size = 18;
n1 = 1.990986606e-7;
a = 1.495978714e8; % kilometers
mu = 3.040357143e-6;

% Select desired Libration point

libind = 2; % 1 for L1 2 for L2 no L3 available

libchange = [-1, 1];

% Find Location of L1 L2 based on roots of polinomial 
L1 = roots([1,-(3-mu),(3-2*mu),-mu,2*mu,-mu]);
L2 = roots([1,(3-mu),(3-2*mu),-mu,-2*mu,-mu]);

% x-position of Lagrange points in rotating frame
% Note xL1 is non dimensional
x_l1 = L1(imag(L1)==0 & real(L1)>0);
x_l2 = L2(imag(L2)==0 & real(L2)>0);

% distance from Sun (mass 1 - primary) to Lagrange points 1, 2
% Non dimensionalized unit
r1_L1 = abs(1-x_l1);
r1_L2 = abs(1+x_l2);


%% richarson

Az = 110000; % kilometers

switch libind
    case 1
        % L1 
        gamma_l = (x_l1);
    case 2
        gamma_l = (x_l2);
end

% dimensionalized quantity stating distance to Libration point
rE = gamma_l*a; 
Az = Az/rE;

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

tau = 0:0.01:2*pi;

tau1 = lamb*w*tau + phi;

deltan = 2 - n;

x = a21*Ax^2 + a22*Az^2 - Ax*cos(tau1) + (a23*Ax^2 - a24*Az^2)*cos(2*tau1) + ...
    (a31*Ax^3 - a32*Ax*Az^2)*cos(3*tau1);

y = k*Ax*sin(tau1) + (b21*Ax^2 - b22*Az^2)*sin(2*tau1) + ...
    (b31*Ax^3 - b32*Ax*Az^2)*sin(3*tau1);

z = deltan*Az*cos(tau1) + deltan*d21*Ax*Az*(cos(2*tau1) - 3) + ...
    deltan*(d32*Az*Ax^2 - d31*Az^3)*cos(3*tau1);

xd = w*lamb*(Ax*sin(tau1) - 2*(a23*Ax^2 - a24*Az^2)*sin(2*tau1) - ...
    3*(a31*Ax^3 - a32*Ax*Az^2)*sin(3*tau1));

yd = w*lamb*(k*Ax*cos(tau1) + 2*(b21*Ax^2 - b22*Az^2)*cos(2*tau1) + ...
    3*(b31*Ax^3 - b32*Ax*Az^2)*cos(3*tau1));

zd = w*lamb*(-deltan*Az*sin(tau1) - 2*deltan*d21*Ax*Az*sin(2*tau1) - ... 
    3*deltan*(d32*Az*Ax^2 - d31*Az^3)*sin(3*tau1));


xlin = x*rE;
ylin = y*rE;
zlin = z*rE;

% figure

figure
% Y Z Axis
subplot(2,2,1)
plot(y, z,'k-')
xlabel('$y$ [$\gamma_L$]' , 'Interpreter', 'Latex', 'FontSize', font_size)
ylabel('$z$ [$\gamma_L$]' , 'Interpreter', 'Latex', 'FontSize', font_size)
set(gca,'fontsize',font_size-4)
grid on
axis equal
% X Z axis
subplot(2,2,2)
plot(x, z,'k-')
xlabel('$x$ [$\gamma_L$]' , 'Interpreter', 'Latex', 'FontSize', font_size)
ylabel('$z$ [$\gamma_L$]' , 'Interpreter', 'Latex', 'FontSize', font_size)
set(gca,'fontsize',font_size-4)
grid on
axis equal

% X Y axis
subplot(2,2,3)
plot(x, y,'k-')
xlabel('$x$ [$\gamma_L$]' , 'Interpreter', 'Latex', 'FontSize', font_size)
ylabel('$y$ [$\gamma_L$]' , 'Interpreter', 'Latex', 'FontSize', font_size)
set(gca,'fontsize',font_size-4)
grid on
axis equal
% 3d
subplot(2,2,4)
plot3(x, y, z,'k-')
xlabel('$x$ [$\gamma_L$]' , 'Interpreter', 'Latex', 'FontSize', font_size)
ylabel('$y$ [$\gamma_L$]' , 'Interpreter', 'Latex', 'FontSize', font_size)
zlabel('$z$ [$\gamma_L$]' , 'Interpreter', 'Latex', 'FontSize', font_size)
set(gca,'fontsize',font_size-4)
grid on
axis equal
view([-41 21])

% non dimensionalized (gamma_l) values for initial conditions
x0 = x(1); y0 = y(1); z0 = z(1); xd0 = xd(1); yd0 = yd(1); zd0 = zd(1);

% vector
x0_l = [x0, y0, z0, xd0, yd0, zd0]';

% initial vector in 1/a1 values (non dimensional
x0_ = (x0_l*rE + [1 - mu + libchange(libind)*gamma_l, 0, 0, 0, 0, 0]'*a)/a;

%% Testing Richardson solution

opts = odeset('RelTol',1e-14,'AbsTol',1e-14);
tspan = [0, pi];

[t,xnumvec] = ode45(@(t,x) eom_3body_rot_fun(t,x,mu), tspan, x0_, opts);

% xnum = (xnumvec(:,1)-(1 - mu + libchange(libind)*gamma_l))/gamma_l;
% ynum = xnumvec(:,2)/gamma_l;
% znum = xnumvec(:,3)/gamma_l;

% Halo in km
xa = xlin + (1-mu+gamma_l)*a;
ya = ylin;
za = zlin;

xnum = xnumvec(:,1)*a;
ynum = xnumvec(:,2)*a;
znum = xnumvec(:,3)*a;

figure
% Y Z Axis
subplot(2,2,1)
plot(ynum, znum,'k-')
hold on
plot(ya,za, 'b--')
xlabel('$y$ [Km]' , 'Interpreter', 'Latex', 'FontSize', font_size)
ylabel('$z$ [Km]' , 'Interpreter', 'Latex', 'FontSize', font_size)
set(gca,'fontsize',font_size-4)
grid on
axis equal
% X Z axis
subplot(2,2,2)
plot(xnum, znum,'k-')
hold on
plot(xa, za,'b--')
xlabel('$x$ [Km]' , 'Interpreter', 'Latex', 'FontSize', font_size)
ylabel('$z$ [Km]' , 'Interpreter', 'Latex', 'FontSize', font_size)
set(gca,'fontsize',font_size-4)
grid on
axis equal

% X Y axis
subplot(2,2,3)
plot(xnum, ynum,'k-')
hold on
plot(xa, ya,'b--')
xlabel('$x$ [Km]' , 'Interpreter', 'Latex', 'FontSize', font_size)
ylabel('$y$ [Km]' , 'Interpreter', 'Latex', 'FontSize', font_size)
set(gca,'fontsize',font_size-4)
grid on
axis equal
% 3d
subplot(2,2,4)
plot3(xnum, ynum, znum,'k-')
hold on
plot3(xa, ya, za,'b--')
xlabel('$x$ [Km]' , 'Interpreter', 'Latex', 'FontSize', font_size)
ylabel('$y$ [Km]' , 'Interpreter', 'Latex', 'FontSize', font_size)
zlabel('$z$ [Km]' , 'Interpreter', 'Latex', 'FontSize', font_size)
set(gca,'fontsize',font_size-4)
grid on
axis equal
view([-41 21])

%% Testing differential corrector
% close all

xstate_init = x0_;

STM_0 = eye(6);
figure

iter = 0;

tol = 1e-12;
ErrTol = 1;
while ErrTol > tol || iter >= 20

xcombined_init = [xstate_init; STM_0(:)];

opts = odeset('RelTol',1e-13,'AbsTol',1e-14, 'Event', @diffcorr_3body_zx_cross_fun);
tspan = [0, 3*1.50623074975922];

[t,xcombinednumvec, te,ye, ie] = ode45(@(t,x) eom_3bodystate_STM_fun(t,x,mu), tspan, xcombined_init, opts);

% Unpack end states
x0 = xcombinednumvec(1, 1:6)';
xf = xcombinednumvec(end, 1:6)';

% Compute dXvec/xt at xf
xdot = eom_3body_rot_fun(0, xf, mu);

phi_tf_t0 = reshape(xcombinednumvec(end, 7:end), 6,6);

% comptue deltax = phi*x0 + dx(tf)/xt * d(T/2);

% % Modification of both dx dz and dyd
% phi_3d = [phi_tf_t0(4,1) - xdot(4)/xf(5)*phi_tf_t0(2,1), phi_tf_t0(4,3) - ...
%     xdot(4)/xf(5)*phi_tf_t0(2,3), phi_tf_t0(4,5) - xdot(4)/xf(5)*phi_tf_t0(2,5);
%     phi_tf_t0(6,1) - xdot(6)/xf(5)*phi_tf_t0(2,1), phi_tf_t0(6,3) - ...
%     xdot(6)/xf(5)*phi_tf_t0(2,3), phi_tf_t0(6,5) - xdot(6)/xf(5)*phi_tf_t0(2,5)];
% 
% delta = pinv(phi_3d)*[-xf(4); -xf(6)];

% Alternative:
%  set dz = 0
phi_3d = [phi_tf_t0(4,1) - xdot(4)/xf(5)*phi_tf_t0(2,1), phi_tf_t0(4,5) - xdot(4)/xf(5)*phi_tf_t0(2,5);
    phi_tf_t0(6,1) - xdot(6)/xf(5)*phi_tf_t0(2,1), phi_tf_t0(6,5) - xdot(6)/xf(5)*phi_tf_t0(2,5)];

% Delta Vector
delta = (phi_3d)\[-xf(4); -xf(6)];

xstate_init(1) = xstate_init(1) + delta(1);
% xstate_init(3) = xstate_init(3) + delta(2);
xstate_init(5) = xstate_init(5) + delta(2);

if 1 

xnum = (xcombinednumvec(:,1) - (1 - mu + libchange(libind)*gamma_l))/gamma_l;
ynum = xcombinednumvec(:,2)/gamma_l;
znum = xcombinednumvec(:,3)/gamma_l;
    
% Y Z Axis
subplot(2,2,1)
plot(ynum, znum,'k-')
hold on
plot(y,z, 'b--')
xlabel('$y$ [$\gamma_L$]' , 'Interpreter', 'Latex', 'FontSize', font_size)
ylabel('$z$ [$\gamma_L$]' , 'Interpreter', 'Latex', 'FontSize', font_size)
set(gca,'fontsize',font_size-4)
grid on
% axis equal
% X Z axis
subplot(2,2,2)
plot(xnum, znum,'k-')
hold on
plot(x, z,'b--')
xlabel('$x$ [$\gamma_L$]' , 'Interpreter', 'Latex', 'FontSize', font_size)
ylabel('$z$ [$\gamma_L$]' , 'Interpreter', 'Latex', 'FontSize', font_size)
set(gca,'fontsize',font_size-4)
grid on
% axis equal

% X Y axis
subplot(2,2,3)
plot(xnum, ynum,'k-')
hold on
plot(x, y,'b--')
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
plot3(x, y, z,'b--')
xlabel('$x$ [$\gamma_L$]' , 'Interpreter', 'Latex', 'FontSize', font_size)
ylabel('$y$ [$\gamma_L$]' , 'Interpreter', 'Latex', 'FontSize', font_size)
zlabel('$z$ [$\gamma_L$]' , 'Interpreter', 'Latex', 'FontSize', font_size)
set(gca,'fontsize',font_size-4)
grid on
% axis equal
view([-41 21])
end
drawnow

iter = iter +1;
ErrTol = sqrt(abs(xf(4))^2 + abs(xf(6))^2);
end
fprintf(['Differential Corrector Succeded, Vel. error ' num2str(ErrTol) '\n'])
%% compute final orbit

opts = odeset('RelTol',1e-14,'AbsTol',1e-14);
tspan = [0, 4.5*pi];

[t,xnumvec] = ode45(@(t,x) eom_3body_rot_fun(t,x,mu), tspan, xstate_init, opts);

xnum = (xnumvec(:,1)-(1 - mu + libchange(libind)*gamma_l))/gamma_l;
ynum = xnumvec(:,2)/gamma_l;
znum = xnumvec(:,3)/gamma_l;

% Halo in km
xa = (xlin + (1-mu+gamma_l)*a)/a;
ya = ylin/a;
za = zlin/a;

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
axis equal
% X Z axis
subplot(2,2,2)
plot(xnum, znum,'k-')
hold on
plot(xa, za,'b--')
xlabel('$x$ [Au]' , 'Interpreter', 'Latex', 'FontSize', font_size)
ylabel('$z$ [Au]' , 'Interpreter', 'Latex', 'FontSize', font_size)
set(gca,'fontsize',font_size-4)
grid on
axis equal

% X Y axis
subplot(2,2,3)
plot(xnum, ynum,'k-')
hold on
plot(xa, ya,'b--')
xlabel('$x$ [Au]' , 'Interpreter', 'Latex', 'FontSize', font_size)
ylabel('$y$ [Au]' , 'Interpreter', 'Latex', 'FontSize', font_size)
set(gca,'fontsize',font_size-4)
grid on
axis equal
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
axis equal
view([-41 21])

geoalt = 42164/a;
moonalt = 384400/a;
geo_orb = [geoalt*cos(t), geoalt*sin(t), 0*cos(t)];
moon_orb = [moonalt*cos(t), moonalt*sin(t), 0*cos(t)];

figure; 
plot3(1, 0, 0, '*k')
hold on
plot3(1+ geo_orb(:,1), geo_orb(:,2), geo_orb(:,3), '--k')
plot3(1+ moon_orb(:,1), moon_orb(:,2), moon_orb(:,3), '--b')
plot3(xa, ya, za,'b--')
plot3(xnum, ynum, znum,'k-')
xlabel('$x$ [Km]' , 'Interpreter', 'Latex', 'FontSize', font_size)
ylabel('$y$ [Km]' , 'Interpreter', 'Latex', 'FontSize', font_size)
zlabel('$z$ [Km]' , 'Interpreter', 'Latex', 'FontSize', font_size)
set(gca,'fontsize',font_size-4)
grid on
view([-16, 5])
leg = legend({'Earth', 'GEO','Moon', 'linear halo', 'nonlinear halo'}, 'Interpreter', 'Latex');
leg.FontSize = font_size;
axis equal

