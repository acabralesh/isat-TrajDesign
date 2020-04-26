%% Test parameter for controller

close all;

% prepare parameters
params_sim;

% define t vector

tvec = 0:SIM.TIME.Ts:(92.4*MINUTES);

lentvec = length(tvec);
%% Sim data out

sim_out.tvec = tvec;
% SC1 data 
sim_out.SC1.x_eci = zeros(lentvec, 2*SIM.UNI.ndim);
% SC2 data
sim_out.SC2.x_eci = zeros(lentvec, 2*SIM.UNI.ndim);
sim_out.SC2.x_lvlh = zeros(lentvec, 2*SIM.UNI.ndim);
sim_out.SC2.x_error = zeros(lentvec, 2*SIM.UNI.ndim);
sim_out.SC2.delta_v = zeros(lentvec, SIM.UNI.ndim);
sim_out.SC2.force_lvlh = zeros(lentvec, SIM.UNI.ndim);

%% Commanded positions 

% Command a natural unstable orbit
SC2.CMD.pos_lvlh = [-2, 0, 0]' * METERS;
SC2.CMD.vel_lvlh = [0, 0, 0]' * METERS/SECONDS;

x_cmd = [SC2.CMD.pos_lvlh; SC2.CMD.vel_lvlh];

%% Initiate sim vectors

x_sc1_eci = [SC1.init_pos_eci; SC1.init_vel_eci];
x_sc2_eci = [SC2.init_pos_eci; SC2.init_vel_eci];

for ii = 1:lentvec
    
    % convert to lvlh frame
    [pos_sc2_lvlh,vel_sc2_lvlh, ~, lvlh_Q_eci] = eci_to_lvlh(SIM.UNI.mu,...
        x_sc1_eci(1:SIM.UNI.ndim), x_sc1_eci(SIM.UNI.ndim+1:2*SIM.UNI.ndim),...
        x_sc2_eci(1:SIM.UNI.ndim), x_sc2_eci(SIM.UNI.ndim+1:2*SIM.UNI.ndim));
    
    % compute error
    
    x_sc2_lvlh = [pos_sc2_lvlh; vel_sc2_lvlh];
    x_error =  x_cmd - x_sc2_lvlh;
    
    delta_v = zeros(SIM.UNI.ndim,1);
            
    force_lvlh = zeros(SIM.UNI.ndim,1);
            
    switch SC2.CTRL.mode
        
        case SC2.CTRL.no_control
            
        case SC2.CTRL.impulsive
            
            u = SC2.CTRL.Klqr*x_error;
            
            delta_v = u;
            
            % need to convert to ECI
            
            
        case SC2.CTRL.zoh
            
            u = SC2.CTRL.Klqr*x_error;
            
            force_lvlh = u;
                  
    end
    
    % change the force into eci frame
    force_eci = lvlh_Q_eci'*force_lvlh;
    
    % update the values in lvlh frame    
    vel_sc2_lvlh = vel_sc2_lvlh + delta_v;
    
    % convert to eci
    [pos_sc2_eci, vel_sc2_eci, eci_Q_lvlh] = lvlh_to_eci(x_sc1_eci(1:SIM.UNI.ndim),...
        x_sc1_eci(SIM.UNI.ndim+1:2*SIM.UNI.ndim),pos_sc2_lvlh,vel_sc2_lvlh);
    
    x_sc2_eci = [pos_sc2_eci; vel_sc2_eci];
    
    % save data
    
    % SC1 data
    sim_out.SC1.x_eci(ii,:) = x_sc1_eci';
    % SC2 data
    sim_out.SC2.x_eci(ii,:) = x_sc2_eci';
    sim_out.SC2.x_lvlh(ii,:) = x_sc2_lvlh';
    sim_out.SC2.x_error(ii,:) = x_error';
    sim_out.SC2.delta_v(ii,:) = delta_v';
    sim_out.SC2.force_lvlh(ii,:) = force_lvlh;
    
    
    % propagate SC1    
    [t,xvec_sc1] = ode45(@(t,x) odekep(t,x,0*force_eci,SIM.UNI.mu,SC1.PHY.mass), ...
        [0 SIM.TIME.Ts], x_sc1_eci, SIM.ODE.opts);
    
    [t,xvec_sc2] = ode45(@(t,x) odekep(t,x,force_eci,SIM.UNI.mu,SC2.PHY.mass), ...
        [0 SIM.TIME.Ts], x_sc2_eci, SIM.ODE.opts);
    
    % rename vectors
    x_sc1_eci = xvec_sc1(end, :)';
    
    x_sc2_eci = xvec_sc2(end, :)';
   
end


%% Plot stuff
font_size = SIM.GRH.font_size;

% plot the states
plot_lvlh_states(sim_out.SC2.x_lvlh,SIM.GRH.font_size)

% Plot the controls
plot_lvlh_control(sim_out.tvec,sim_out.SC2, SC2, font_size)

% Plot state error
plot_lvlh_error(sim_out.tvec,sim_out.SC2.x_error, SC2, font_size)

