%% Main drone parameters and filters

%% Simulation parameters
Ts    = 1/1000;     % Sample time PPRZ controller, reference trajectories 
t_simu           = 25;
sample_t_simu    = 1e-3;      % Sample time simulation

load("models/slexAircraftPitchControlData.mat");



% 
% [A,B,C,D] = linmod("AC_z");
% 
% AC_z = ss(A,B,C,D);
% 
% 
% pzmap(AC_z)