% Clean up the Manifolds script
% load('EM_results_2_2.mat');
load('results_2_18.mat');

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

% mu = M_moon/(M_moon + M_earth);

if 0

% EML2
EML1_manifolds.mu = mu;
EML1_manifolds.a = a;
EML1_manifolds.description = 'Manifolds from EML1 with amplitude of 30000. All value are unitless';

% select first 45 orbits
for ii = 1:40
    orb_ind = results.orb_sort_ind(ii);
    
    xnumvec = results.sol{orb_ind}.xnumvec;
    tnumvec = results.sol{orb_ind}.t;
    
    tmin_id = results.minind(orb_ind);
    
    xnumvec = xnumvec(1:tmin_id, :);
    tnumvec = tnumvec(1:tmin_id);
    
    EML1_manifolds.sol{ii}.xnumvec = xnumvec;
    EML1_manifolds.sol{ii}.tnumvec = tnumvec;
    
end

elseif 1
    
    
    % EML2
SEL2_manifolds.mu = mu;
SEL2_manifolds.a = a;
SEL2_manifolds.description = 'Manifolds from SEL2 with amplitude of 110000. All value are unitless';

% select first 45 orbits
for ii = 1:35
    orb_ind = results.orb_sort_ind(ii);
    
    xnumvec = results.sol{orb_ind}.xnumvec;
    tnumvec = results.sol{orb_ind}.t;
    
    tmin_id = results.minind(orb_ind);
    
    xnumvec = xnumvec(1:tmin_id, :);
    tnumvec = tnumvec(1:tmin_id);
    
    SEL2_manifolds.sol{ii}.xnumvec = xnumvec;
    SEL2_manifolds.sol{ii}.tnumvec = tnumvec;
    
end
end

