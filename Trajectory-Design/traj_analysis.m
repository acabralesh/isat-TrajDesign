% Close all previous figures and clear terminal and workspace
close all
clc
clear all

% Add the trajectory design directory
addpath(genpath('.'));  % main directory
addpath(genpath('../Telescope-Parameters/Units'));

%% Launch Parameters
% parameters for a SINGLE PRIORITY PHASE
Model.LAN.t_step = 1; % number of days in a timestep, days
Model.LAN.t_start = [2021 01 01]; % [year month day] e.g. [2021 6 18]
Model.LAN.t_end = [2021 1 31]; % [year month day] e.g. [2021 8 18]

%% Determine the parameters for the optimal orbit for all launches
% Launch Location
input_struct.launch_lat = 28*DEGREES; 
input_struct.launch_long = 20*DEGREES;
% Atmospheric Drag
input_struct.atmosphere_drag = 1.5*KILOMETERS/SECONDS; % 
% Print stuff
input_struct.verbose = 1;

%%%%%---------------------------------------------------------
%  Assembly location parameter
%      ---KEY---
%          1: Earth-Moon L1 (Lunar Gateway)
%          2: Sun-Earth L2
%%%%%---------------------------------------------------------
locationOptions = ["EML1","SEL2"];
assembly_loc = locationOptions(1);

delta_v = [];
t_launch = {};
t_duration = [];

% Find launch time options at t_step intervals
for t_sample = datetime(Model.LAN.t_start):caldays(Model.LAN.t_step):datetime(Model.LAN.t_end)
    [delta_v(end+1), t_launch{end+1}, t_duration(end+1), ~] = traj_lookup(t_sample, t_sample, assembly_loc, input_struct);
end

% Plot delta-v for the different trajectories
figure
hold on
days = 1:length(delta_v);
subplot(2,1,1);
plot(days, delta_v)
xlim([0 length(delta_v)])
xlabel('Days','FontSize',12)
ylabel('Delta-V','FontSize',12)
subplot(2,1,2);
plot(days, t_duration)
xlim([0 length(delta_v)])
xlabel('Days','FontSize',12)
ylabel('Duration','FontSize',12)
