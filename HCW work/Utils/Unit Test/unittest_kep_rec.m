% Unit Test for Kepler/Rectangular Transforms

Apophis_pos = [0.2991634567614704,
               0.869437579145023,
               -0.039177781942555084];
Apophis_vel = [-0.015541893965049815,
             0.008948110077735555,
            -0.0008458270588252229];

keplerconstant = 0.01720209895^2;

[a, eccentricity, inclination, longnode, argperi, meananom] =...
 rec_to_kepler(keplerconstant, Apophis_pos, Apophis_vel);


[pos,vel] = kepler_to_rec(keplerconstant,a, eccentricity, inclination, longnode, argperi, meananom);


Apophis_pos-pos

Apophis_vel-vel