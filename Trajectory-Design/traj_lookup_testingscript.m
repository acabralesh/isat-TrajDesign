%% Testing traj_lookup function

%% Input flag parameters currently inplemented
% Launch Location
input_struct.launch_lat = 28*DEGREES; 
input_struct.launch_long = 20*DEGREES;
% Atmospheric Drag
input_struct.atmosphere_drag = 1.5*KILOMETERS/SECONDS; % 
% Print stuff
input_struct.verbose = 1;



%% Run optimization look
t_start = [2021 6 18];
t_end = [2021 6 19];
target_loc = 'SEL2';
[delta_v, t_launch, t_duration, out_struct] = traj_lookup(t_start, t_end, target_loc, input_struct);

delta_v
datetime(t_launch)