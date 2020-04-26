% Test the generation of the rotation matrix


target_loc = 'SEL2';

JD = juliandate(2000, 6,15);

switch target_loc
    
    case 'SEL2'
        
        [pos, vel] = planetEphemeris(JD, 'Sun', 'Earth');        
        
        % convert into column vector
        pos = pos(:);
        vel = vel(:);
        
        i_hat = pos/norm(pos);
        k_hat = cross(pos, vel)/norm(cross(pos, vel));
        j_hat = cross(k_hat, i_hat);
        
        cer_C_ber = [i_hat, j_hat, k_hat];
        
    case 'EML1'
        [pos, vel] = planetEphemeris(JD, 'Earth', 'Moon');        
        
        % convert into column vector
        pos = pos(:);
        vel = vel(:); 
        
        i_hat = pos/norm(pos);
        k_hat = cross(pos, vel)/norm(cross(pos, vel));
        j_hat = cross(k_hat, i_hat);
        
        cer_C_ber = [i_hat, j_hat, k_hat];
        
    otherwise
        error('Currently only SEL2 and EML1 system are supported');
end


load('results.mat')

xnumvec = results.sol{12}.xnumvec;
tnumvec = results.sol{12}.t;


% Generate the position in ECEF
lat = 28*DEGREES;
long = 20*DEGREES;
launch_pos_ecef =lla2ecef( [lat/DEGREES, long/DEGREES, 0]);



earth_radius = norm(launch_pos_ecef);

[r_ECI, v_ECI, a_ECI] = ECEFtoECI(JD,launch_pos_ecef',[0 0 0]',[0,0,0]');


figure;
hold on
plot_earth(earth_radius, 0)
plot3(r_ECI(1),r_ECI(2),r_ECI(3), '*k')
quiver3(r_ECI(1),r_ECI(2),r_ECI(3),1e3*v_ECI(1), 1e3*v_ECI(2), 1e3*v_ECI(3),'k-', 'linewidth', 2)

