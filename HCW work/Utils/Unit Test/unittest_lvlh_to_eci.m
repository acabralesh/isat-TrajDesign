% unittest_lvlh_to_eci

RA = [-266.77,  3865.8, 5426.2];     % [km]
VA = [-6.4836, -3.6198, 2.4156];     % [km/s]
RB = [-5890.7, -2979.8, 1792.2];     % [km]
VB = [0.93583, -5.2403, -5.5009];    % [km/s]


% generate one way ECI -> LVLH
% Earth gravitational parameter
mu = 398600;                        % [km^3/s^2]

[SC2_rel_pos_LVLH,SC2_rel_vel_LVLH, SC1_rel_a_LVLH, lvlh_Q_eci] = eci_to_lvlh(mu, ...
    RA', VA', RB', VB');

% other way LVLH -> ECI

[RB2,VB2] = lvlh_to_eci(RA', VA', SC2_rel_pos_LVLH,SC2_rel_vel_LVLH);

% difference
RB'-RB2

VB'-VB2