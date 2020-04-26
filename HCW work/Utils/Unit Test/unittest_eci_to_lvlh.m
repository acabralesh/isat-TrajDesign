% lvlhtestscript
% 
% % kepler constant
% mu = 3.986004418*10^(14)*METERS^3/SECONDS^2;
% 


RA = [-266.77,  3865.8, 5426.2];     % [km]
VA = [-6.4836, -3.6198, 2.4156];     % [km/s]
RB = [-5890.7, -2979.8, 1792.2];     % [km]
VB = [0.93583, -5.2403, -5.5009];    % [km/s]

% Earth gravitational parameter
mu = 398600;                        % [km^3/s^2]

[SC2_rel_pos_LVLH,SC2_rel_vel_LVLH, SC1_rel_a_LVLH, lvlh_Q_eci] = eci_to_lvlh(mu, ...
    RA', VA', RB', VB')
