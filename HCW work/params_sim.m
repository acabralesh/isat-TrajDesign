%% ----------------------Spacecraft Properties File----------------------
addpath('Utils')
addpath('../Units')

%% Graphics Properties
SIM.GRH.font_size = 18;

%% Universal Properties
SIM.UNI.G = 6.67408e-11 * METERS^3 * KILOGRAMS^-1 * SECONDS^-2;
SIM.UNI.M_sun = 1.9884754e30 * KILOGRAMS;
SIM.UNI.M_earth = 5.972358e24 * KILOGRAMS;
SIM.UNI.M_moon = 7.34767309e22 * KILOGRAMS;
% Gravitational parameter
SIM.UNI.mu = SIM.UNI.G*SIM.UNI.M_earth;

SIM.UNI.earthradius = 6371*KILOMETERS;

SIM.UNI.ndim = 3;

SIM.TIME.Ts = 10 * SECONDS;

%% ODE properties
SIM.ODE.opts = odeset('RelTol',1e-10,'AbsTol',1e-10);
%% SC1 Orbital Parameters
SC1.ORB.a = SIM.UNI.earthradius + 400*KILOMETERS;
SC1.ORB.e = 0;
SC1.ORB.inc = 0 * DEGREES; 
SC1.ORB.longnode = 0 * DEGREES;
SC1.ORB.argperi = 0 * DEGREES;
SC1.ORB.meananom = 0 * DEGREES;

% mean orbital rate
SC1.ORB.n = sqrt(SIM.UNI.mu/SC1.ORB.a^3);

% Convert to ECI frame
[SC1.init_pos_eci, SC1.init_vel_eci] = kepler_to_rec(SIM.UNI.mu, SC1.ORB.a, ...
    SC1.ORB.e, SC1.ORB.inc, SC1.ORB.longnode, SC1.ORB.argperi, SC1.ORB.meananom);

%% SC1 Physical Properties
SC1.PHY.mass = 1 * KILOGRAMS; 

%% SC2 Physical Properties
SC2.PHY.mass = 1 * KILOGRAMS;

%% SC2 Initial Position
SC2.init_pos_lvlh = [0, -2.5, 0]' * METERS;
SC2.init_vel_lvlh = [0, 0, 0]' * METERS/SECONDS;

% convert to eci frame
[SC2.init_pos_eci,SC2.init_vel_eci, SC2.init_eci_Q_lvlh] = ...
    lvlh_to_eci(SC1.init_pos_eci,  SC1.init_vel_eci,...
    SC2.init_pos_lvlh,SC2.init_vel_lvlh);

%% SC2 Controller Options

% Control sampling time
SC2.CTRL.Ts = SIM.TIME.Ts; %10 * SECONDS; 

SC2.CTRL.no_control = 0;
SC2.CTRL.impulsive = 1;
SC2.CTRL.zoh = 2;

% select the controller mode
SC2.CTRL.mode = SC2.CTRL.no_control;
%% SC2 Controller gains

[SC2.CTRL.Ac, SC2.CTRL.Bc] = conmat_hcw(SC1.ORB.n, SIM.UNI.ndim);

SC2.CTRL.gamma = 1e-3;

SC2.CTRL.Q = SC2.CTRL.gamma*eye(2*SIM.UNI.ndim);

SC2.CTRL.R = eye(SIM.UNI.ndim);

switch SC2.CTRL.mode
    
    case SC2.CTRL.no_control
        % set the correct LQR gains
        SC2.CTRL.Klqr = zeros(SIM.UNI.ndim, 2*SIM.UNI.ndim);
        
    case SC2.CTRL.impulsive
        
        SC2.CTRL.Phi = expm(SC2.CTRL.Ac*SC2.CTRL.Ts) + 0*stm_HCW(SC1.ORB.n, SC2.CTRL.Ts);
        
        SC2.CTRL.Ad = SC2.CTRL.Phi;
        
        % for impulsive system B = [Phi_rv; Phi_vv];
        SC2.CTRL.Bd = SC2.CTRL.Phi(:, SIM.UNI.ndim+1:end);
        
        SC2.CTRL.Klqr = dlqr(SC2.CTRL.Ad, SC2.CTRL.Bd, SC2.CTRL.Q, SC2.CTRL.R);
        
    case SC2.CTRL.zoh
        
        % zero order hold use lqrd
        
        SC2.CTRL.Klqr = lqrd(SC2.CTRL.Ac, SC2.CTRL.Bc, ...
            SC2.CTRL.Q, SC2.CTRL.R, SC2.CTRL.Ts);
end


