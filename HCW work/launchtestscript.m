% Launch orbit test
% close all

addpath('Utils')
addpath('../Units')
mu = 3.986004418*10^(14)*METERS^3/SECONDS^2;

earth_radius = 6351*KILOMETERS;

opts = odeset('RelTol',1e-10,'AbsTol',1e-10);


r_orb = earth_radius + 500*KILOMETERS;

v_earth = 0*491*METERS/SECONDS;

dv_launch = sqrt(mu/earth_radius)*(sqrt(2*r_orb/(r_orb + earth_radius)));

xinit = [earth_radius; 0; 0; 0*1000; (dv_launch+v_earth); 0]

T = 2*pi*sqrt(r_orb^3/mu);

t = [0 2*T];


[t,xvec_sc1] = ode45(@(t,x) odekep(t,x,0,mu,1), ...
        t, xinit, opts);
    
figure;
plot(xvec_sc1(:,1), xvec_sc1(:,2))
hold on
plot(xvec_sc1(end,1), xvec_sc1(end,2), '*')
grid on
axis square
