function [v1_lamvec,v2_lamvec, T_lamb] = lambert_mina(mu,x1_state, x2_state)
%[v1_lamvec,v2_lamvec, T_lamb] = LAMBERT_MINA(mu,x1_state,
%x2_state)
%
% Purpose: Computes the delta-v to solve the minimum ellipse lambert's
% problem
%
% Inputs:
% mu        [1x1] Gravitational parameter (mass ratio parameter) M2/(M1+M2)
%           given that M1 > M2
% x1_state  [6x1] State of initial position and velocityfor the boundary
%           value
% x2_state  [6x1] State of final position and velocity for the boundary
%           value
%
% Outputs: 
% v1_lamvec [3x1] Initial velocity for the transfer trajectory
% v2_lamvec [3x1] Final velocity for the transfer trajectory
% T_lamb    [1x1] Transfer time
%
% Co-dependencies
%  None
%
% Source:
%  [1] Prussing, John E.    Orbital mechanics / John E. Prussing, Professor
%  Emeritus of Aerospace Engineering, University of Illinois at
%  Urbana-Champaign, Bruce A. Conway, Professor of Aerospace Engineering,
%  University of Illinois at Urbana-Champaign.  New York : Oxford
%  University Press, [2013]
%
% Author: Alejandro Cabrales H


r1_vec = x1_state(1:3);
v1_vec = x1_state(4:6);
r2_vec = x2_state(1:3);
v2_vec = x2_state(4:6);

r1 = norm(r1_vec);
r2 = norm(r2_vec);

n_vec = cross(r1_vec, v1_vec);

e_n = n_vec/norm(n_vec);

e_d = dot(r1_vec/r1, cross(r2_vec/r2,e_n));

theta = acos(dot(r1_vec,r2_vec)/(r1*r2));

if sign(e_d) >= 0
%     theta = theta;
else
    theta = 2*pi-theta;
end
% chord between two points
c_vec = r2_vec-r1_vec;

c = norm(c_vec);

% perimeter of space triangle
s = 1/2*(r1 + r2 + c);

u1 = r1_vec/norm(r1);

u2 = r2_vec/norm(r2);

uc = c_vec/c;

% minimum semi major axis
am = s/2; 

% define lamberts alpha betas
alpha_m = pi;

beta_m = 2*asin(((s-c)/s)^(1/2));

% time associated with transfer
tm = 1/sqrt(mu) * (s^3/8)^(1/2)*(alpha_m - beta_m + sin(beta_m));

% ---- conditions on alpha beta ------

% select a
a = am;

% solve for alpha_0, beta_0

alpha_0 = 2*asin((s/(2*a))^0.5);

beta_0 = 2*asin(((s-c)/(2*a))^0.5);

% make the correct beta
if theta < pi
    
    beta = beta_0;
    
else
    beta = -beta_0;
    
end

% ******** deal with the time aspect later on ********
alpha = alpha_0;

A = sqrt(mu/(4*a))*cot(alpha/2);

B = sqrt(mu/(4*a))*cot(beta/2);

v1_lamvec = (B + A)*uc + (B - A)*u1;

v2_lamvec = (B + A)*uc - (B - A)*u2;

% initial position for propagator
% xinit_lamb = [r1_vec; v1_lamvec];

% transfer time (time of flight)
T_lamb = a^(3/2)/sqrt(mu) * ( alpha - beta - (sin(alpha) - sin(beta)));


end

