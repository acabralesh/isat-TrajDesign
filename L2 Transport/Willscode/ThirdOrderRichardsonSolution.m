%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
%   3rd Order Richardson Halo Orbit Approximation
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

format long;
clear; clc; close all;

global params
params.tol = 1e-10;

%% Generate initial conditions from 3rd Order Richardson Solution

atL2 = 1;
atL1 = 0;

n1 = 1.990986606e-7;
a1 = 1.495978714e8;
mu = 3.040357143e-6;
params.Az = 4e4;
params.mu = mu;

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

%%% distance from Jupiter (mass 2 - secondary) to Lagrange points 1, 2, and 3
r2_L1 = abs(xL1);
r2_L2 = abs(xL2);
r2_L3 = abs(1+xL3);

% Energy cases
E1 = -0.5*((1-mu)*r1_L1^2+mu*r2_L1^2)-(1-mu)/r1_L1-mu/r2_L1;
E2 = -0.5*((1-mu)*r1_L2^2+mu*r2_L2^2)-(1-mu)/r1_L2-mu/r2_L2;
E3 = -0.5*((1-mu)*r1_L3^2+mu*r2_L3^2)-(1-mu)/r1_L3-mu/r2_L3;
E4 = -3/2;
E5 = -3/2;

if(atL2)
    gamma = r2_L2;
    rE = r2_L2*a1;
    c2 = (1/gamma^3)*((-1)^2*mu+(-1)^2*(1-mu)*gamma^3/(1+gamma)^3);   %good
    c3 = (1/gamma^3)*((-1)^3*mu+(-1)^3*(1-mu)*gamma^4/(1+gamma)^4);   %good
    c4 = (1/gamma^3)*((-1)^4*mu+(-1)^4*(1-mu)*gamma^5/(1+gamma)^5);   %good
elseif(atL1)
    gamma = r2_L1;
    rE = r2_L1*a1;
    c2 = (1/gamma^3)*((1)^2*mu+(-1)^2*(1-mu)*gamma^3/(1-gamma)^3);   %good
    c3 = (1/gamma^3)*((1)^3*mu+(-1)^3*(1-mu)*gamma^4/(1-gamma)^4);   %good
    c4 = (1/gamma^3)*((1)^4*mu+(-1)^4*(1-mu)*gamma^5/(1-gamma)^5);   %good
end

params.Az = params.Az/rE;

lambda = roots([1,0,(c2-2),0,-(c2-1)*(1+2*c2)]);                  
lambda = lambda(find(imag(lambda)==0 & real(lambda)>0));          %good

wp = sqrt((2-c2+sqrt(9*c2^2-8*c2))/2);                            %good
wv=sqrt(c2);                                                      %good
k = (wp^2+1+2*c2)/(2*wp);                                         %good
delta = wp^2-c2;                                                  %good

d1 = (3*lambda^2/k)*(k*(6*lambda^2-1)-2*lambda);                  %good
d2 = (8*lambda^2/k)*(k*(11*lambda^2-1)-2*lambda);                 %good

a21 = 3*c3*(k^2-2)/(4*(1+2*c2));                                  %good
a22 = 3*c3/(4*(1+2*c2));                                          %good
a23 = (-3*c3*lambda/(4*k*d1))*(3*k^3*lambda-6*k*(k-lambda)+4);    %good
a24 = (-3*c3*lambda/(4*k*d1))*(2+3*k*lambda);                     %good

d21 = -c3/(2*lambda^2);                                           %good
d31 = (3/(64*lambda^2))*(4*c3*a24+c4);                            %good
d32 = (3/(64*lambda^2))*(4*c3*(a23-d21)+c4*(4+k^2));              %good

b21 = (-3*c3*lambda/(2*d1))*(3*k*lambda-4);                       %good
b22 = (3*c3*lambda)/d1;                                           %good
b31 = (3/(8*d2))*(8*lambda*(3*c3*(k*b21-2*a23)-c4*(2+3*k^2))+(9*lambda^2+1+2*c2)*(4*c3*(k*a23-b21)+k*c4*(4+k^2))); %good
b32 = (9*lambda/d2)*(c3*(k*b22+d21-2*a24)-c4) + (3*(9*lambda^2+1+2*c2)/(8*d2))*(4*c3*(k*a24-b22)+k*c4); %good


a31 = (-9*lambda/(4*d2))*(4*c3*(k*a23-b21)+k*c4*(4+k^2)) + ((9*lambda^2+1-c2)/(2*d2))*(3*c3*(2*a23-k*b21)+c4*(2+3*k^2)); %good
a32 = (-9*lambda/(4*d2))*(4*c3*(k*a24-b22)+k*c4)-(3/(2*d2))*(9*lambda^2+1-c2)*(c3*(k*b22+d21-2*a24)-c4);


s1 = (1/(2*lambda*(lambda*(1+k^2)-2*k)))*(1.5*c3*(2*a21*(k^2-2)-a23*(k^2+2)-2*k*b21)-(3/8)*c4*(3*k^4-8*k^2+8)); %good
s2 = (1/(2*lambda*(lambda*(1+k^2)-2*k)))*(1.5*c3*(2*a22*(k^2-2)+a24*(k^2+2)+2*k*b22+5*d21)+(3/8)*c4*(12-k^2));  %good
l1 = -1.5*c3*(2*a21+a23+5*d21)-(3/8)*c4*(12-k^2)+2*lambda^2*s1;   %good
l2 = 1.5*c3*(a24-2*a22)+(9/8)*c4+2*lambda^2*s2;                   %good

Ax = sqrt((-l2*params.Az^2-delta)/l1);
Ax_km = Ax*rE;

w = 1 + s1*Ax^2 + s2*params.Az^2;

T = (2*pi/(lambda*w*n1))/(3600*24);

delta_n = 1; % for northern halo, enter -1 for southern halo
t = 0:0.01:2*pi;
phi = 0;
tau1 = lambda*(w*t)+phi;

xlin = (a21*Ax^2 + a22*params.Az^2 - Ax*cos(tau1) + (a23*Ax^2 - a24*params.Az^2)*cos(2*tau1) + (a31*Ax^3-a32*Ax*params.Az^2)*cos(3*tau1))*rE;
ylin = (k*Ax*sin(tau1) + (b21*Ax^2-b22*params.Az^2)*sin(2*tau1) + (b31*Ax^3 - b32*Ax*params.Az^2)*sin(3*tau1))*rE;
zlin = (delta_n*params.Az*cos(tau1) + delta_n*d21*Ax*params.Az*(cos(2*tau1)-3)+delta_n*(d32*params.Az*Ax^2-d31*params.Az^3)*cos(3*tau1))*rE;

plot3(xlin,ylin,zlin); grid on;

x0 = a21*Ax^2 + a22*(params.Az)^2 - Ax*cos(phi) + (a23*Ax^2 - a24*params.Az^2)*cos(2*phi) + (a31*Ax^3-a32*Ax*params.Az^2)*cos(3*phi);
y0 = k*Ax*sin(phi) + (b21*Ax^2-b22*params.Az^2)*sin(2*phi)+(b31*Ax^3-b32*Ax*params.Az^2)*sin(3*phi);
z0 = delta_n*params.Az*cos(phi) + delta_n*d21*Ax*params.Az*(cos(phi)-3)+delta_n*(d32*params.Az*Ax^2-d31*params.Az^3)*cos(3*phi);
xdot0 = lambda*w*(Ax*sin(phi)-2*(a23*Ax^2-a24*params.Az^2)*sin(2*phi)-3*(a31*Ax^3-a32*Ax*params.Az^2)*sin(3*phi));
ydot0 = lambda*w*(k*Ax*cos(phi) + 2*(b21*Ax^2-b22*params.Az^2)*cos(2*phi) + 3*(b31*Ax^3-b32*Ax*params.Az^2)*cos(3*phi));
zdot0 = -lambda*w*delta_n*(params.Az*sin(phi)+2*d21*Ax*params.Az*sin(2*phi)+3*(d32*params.Az*Ax^2-d31*params.Az^3)*sin(3*phi));

x0_l = [x0;y0;z0;xdot0;ydot0;zdot0];

x0_ = (x0_l*rE+[1-mu+gamma;0;0;0;0;0].*a1)/a1

%% changed alex


opts = odeset('RelTol',1e-14,'AbsTol',1e-14);
[t,z] = ode45('threeDimensional3BodyOrbitRotatingFrame',[0,pi],x0_,opts);
z=z*a1;
figure(1)
plot3(z(:,1),z(:,2),z(:,3)); grid on;
axis equal

xlin = xlin + (1-mu+gamma)*a1;
ylin = ylin;
zlin = zlin;

figure(2)
subplot(2,2,1);
plot(z(:,2),z(:,3), 'linewidth', 1)
hold on
plot(ylin, zlin)
xlabel('y'); ylabel('z')
axis equal
grid on;
subplot(2,2,2);
plot(z(:,1),z(:,3), 'linewidth', 1)
hold on
plot(xlin, zlin)
xlabel('x'); ylabel('z')
axis equal
grid on;
subplot(2,2,3);
plot(z(:,1),z(:,2), 'linewidth', 1)
hold on
plot(xlin, ylin)
xlabel('x'); ylabel('y')
axis equal
grid on;

figure(3)
plot3(z(:,1),z(:,2),z(:,3));
radiusSun = 6.95508e5;
radiusEarth = 6.371e3; 
hold on;
[xs, ys, zs] = sphere;
%surf(xs*radiusSun - mu, ys*radiusSun, zs*radiusSun);
hold on;
surf(xs*radiusEarth + 1 - mu, ys*radiusEarth, zs*radiusEarth);
grid on;






%% Use differential corrections algorithm to obtain nonlinear initial conditions

% Use initial conditions generated from the analytical solution
%   of the linearized system as the initial guess

M0 = [x0_ eye(6)];
[M0_nonlinear] = differentialCorrectionsAlg3D(M0,params.tol);
x0_nonlinear = M0_nonlinear(:,1);

%% Integrate

opts = odeset('RelTol',1e-14,'AbsTol',1e-14);
[t,z] = ode45('threeDimensional3BodyOrbitRotatingFrame',[0,2*pi],x0_nonlinear,opts);
z=z*a1;
figure(1)
plot3(z(:,1),z(:,2),z(:,3)); grid on;
axis equal

figure(2)
subplot(2,2,1);
plot(z(:,2),z(:,3))
axis equal
grid on;
subplot(2,2,2);
plot(z(:,1),z(:,3))
axis equal
grid on;
subplot(2,2,3);
plot(z(:,1),z(:,2))
axis equal
grid on;

figure(3)
plot3(z(:,1),z(:,2),z(:,3));
radiusSun = 6.95508e5;
radiusEarth = 6.371e3; 
hold on;
[xs, ys, zs] = sphere;
%surf(xs*radiusSun - mu, ys*radiusSun, zs*radiusSun);
hold on;
surf(xs*radiusEarth + 1 - mu, ys*radiusEarth, zs*radiusEarth);
grid on;






