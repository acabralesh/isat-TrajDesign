%% Converts BER to ECI

% my values
G = 6.67408e-11 * METERS^3 * KILOGRAMS^-1 * SECONDS^-2;
M_sun = 1.9884754e30 * KILOGRAMS;
M_earth = 5.972358e24 * KILOGRAMS;
M_moon = 7.34767309e22 * KILOGRAMS;

% Sun - Earth
a = 1.4959787e8 * KILOMETERS;
% Earth - Moon
% a = 384400 * KILOMETERS;

% Sun - Earth
mu = M_earth/(M_sun + M_earth);

theta_year = 0;
theta_hour = 0;

% ---- inputs -----
omega_n_ber = [0; 0; 1]; % Angular velocity vector

ecliptic_inc = -23.5*DEGREES;

% initial rotation matrices
bei_R_ber = @(t) [cos(theta_year + t) -sin(theta_year + t) 0;
    sin(theta_year + t) cos(theta_year + t) 0;
    0 0 1];

eci_R_ecef = @(t) [cos(theta_hour + t) -sin(theta_hour + t) 0;
    sin(theta_hour + t) cos(theta_hour + t) 0;
    0 0 1];


cei_R_eci = [1 0 0;
    0 cos(ecliptic_inc) -sin(ecliptic_inc);
    0 sin(ecliptic_inc) cos(ecliptic_inc)];

cei_R_bei = eye(3);
figure
for idclose = 1:length(results.orb_sort_ind)
    
    solid = results.orb_sort_ind(idclose);
    % obtain the minimum results
    
    xnumvec_ber = results.sol{solid}.xnumvec';
    
    t_min_ind = results.minind(solid);
    
    % obtain the minimum time
    tmin = results.sol{solid}.t(t_min_ind);
    
    % obtain the manifold from start to finish
    xnumvec_ber = xnumvec_ber(:, t_min_ind:-1:1);
    time = results.sol{solid}.t(t_min_ind:-1:1) - tmin;
    
    xnumvec_bei = zeros(size(xnumvec_ber));
    xnumvec_cei = zeros(size(xnumvec_ber));
    xnumvec_eci = zeros(size(xnumvec_ber));
    
    
    for ii = 1:length(xnumvec_ber)
        
        % obtain time
        theta = time(ii);
        
        % Convert BER to BEI
        [sc_pos_bei,sc_vel_bei] = ber_to_bei(bei_R_ber(theta), omega_n_ber,xnumvec_ber(1:3,ii), xnumvec_ber(4:6,ii));
        
        xnumvec_bei(:, ii) = [sc_pos_bei; sc_vel_bei];
        
        % Convert BER to CEI
        r_earth_bei = bei_R_ber(theta)*[1-mu;0;0];
        
        [sc_pos_cei,sc_vel_cei] = bei_to_cei(cei_R_bei, omega_n_ber,r_earth_bei, sc_pos_bei, sc_vel_bei);
        
        xnumvec_cei(:, ii) = [sc_pos_cei;sc_vel_cei];
        
%         sc_pos_eci = cei_R_eci'*sc_pos_cei;
%         sc_vel_eci = cei_R_eci'*sc_vel_cei;
        
%         xnumvec_eci(:, ii) = [sc_pos_eci;sc_vel_eci];
        
    end
    
    results.sol{solid}.xnumvec_cei = xnumvec_cei;
    plot3(xnumvec_cei(1,:),xnumvec_cei(2,:),xnumvec_cei(3,:));
    hold on
    
end