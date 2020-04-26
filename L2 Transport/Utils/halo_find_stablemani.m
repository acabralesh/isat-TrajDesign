function sol = halo_find_stablemani(mu, t_stm_span,...
    xstate_init, STM_0, d, t_mani_span, ode_opts)
%HALO_3BODYROT_RICHARDSON_INIT(mu, x_l, libind, Az, phi, n) %
% Purpose: 
%
% Inputs:
% - mu              [1x1] gravitational parameters
% - t_stm_span      [1xn] The time vector that obtains the state transition
%                   matrix. First velue must be equal to 0, it will find
%                   the manifold on the n-1 points given by t_stm_span
% - xstate_init     [6x1] Initial state of the manifold solution given by
%                   richardson solution with Differnetial correction
% - STM_0           [6x6] Initial condition for the state transition matrix
% - d               [1x1] Disturbance distance on the stable eigenvector
%                   direction
% - tmanispan       [1xnt] Time vector that the manifold will be propagated
% - ode_opts        Options for the ode 45 system
%
% Outputs: 
% - sol             Structure with solutions
%   - t             [1xnt] time vector that is the same as t_mani_span
%   - xnumvec       [n-1] vector of state solutions for each manifold 
%
% Co-dependencies
%
% Source:
%  [1] Richardson, D.L. Celestial Mechanics (1980) 22: 241. "Analytic
%  Construction of Periodic Orbits about the Col;inear Points."
%  https://doi.org/10.1007/BF01229511
%  [2] Barde, B. "Application of dynamical systems theory in mission design
%  and conceptual development for libration point missions." MS Thesis.
%  Purdue University, 2000. Chapter 3.
% 
% Author: Alejandro Cabrales H

% solution storage

sol = {};


% set up the combined initial state
xcombined_init = [xstate_init(:); STM_0(:)];

% Compute the orbit plus the state transition matrix
[t,xcombinednumvec] = ode45(@(t,x) eom_3bodystate_STM_fun(t,x,mu), t_stm_span, xcombined_init, ode_opts);


% for each manifold, skip the first one since correspond to 0 0
for ii = 2:size(xcombinednumvec,1)
    
    % obtain the vector
    x_i = xcombinednumvec(ii, 1:6)';
    
    % reshape to obtain phi
    phi = reshape(xcombinednumvec(ii, 7:end), 6,6);
    
    % obtain eigenvectors and associated eigenvalues
    [V,D] = eig(phi);
    Ddiag = diag(D);
    
    % find the unstable (maximum) eigen values (must be imaginary)
    lamb_unst = max(Ddiag(imag(Ddiag) == 0));
    
    % stable eigenvalue
    [lamb_stab, ind_stab] = min(Ddiag - 1/lamb_unst);
     
    eigvec_stab = V(:, ind_stab);
    
    % X_ws obtain the new initial condition
    x_ws = x_i + d*eigvec_stab;

    [t,xnumvec] = ode45(@(t,x) eom_3body_rot_fun(t,x,mu), t_mani_span, x_ws, ode_opts);
    
    sol{ii-1}.xnumvec = xnumvec;
    sol{ii-1}.t = t;
    
end



end

