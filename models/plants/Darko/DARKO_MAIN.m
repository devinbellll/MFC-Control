%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%        DARKO Flight Simulator         %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%          Model : Phi-Theory           %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%    Controller : Model-Free Control    %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;
clc;

%% DarkO initial points for simulation
compute_drone_initial_points;

%% Simulation parameters
t_simu           = 180;
Ts               = 1/500; %512;     % Sample time PPRZ controller, reference trajectories 
sample_t_simu    = 1e-3;            % Sample time simulation
use_control_sat  = 0;
slope_trans      = -0.75;
v_b_b            = [0;0;0];

%% Delta parameters analysis
delta_mass      = 0;
delta_inertia   = 0;
delta_phi_coef  = 0.0;  % Zero for the first analysis IJMAV2019
time_var_param  = 13500;

%% Max delay analysis
delay_max_angles    = 1; % 20 (separated)
delay_max_speeds    = 1; % 65 (separated) 
delay_max_positions = 1; % 90 (separated)

%% Coordinated circle turn
r = 5;
xc = 0;
yc = 0;
time_circle_traj_start = 60;

%% Chose "trim" point i,j
% i = 1  => 90    degrees theta (Hover Flight)
i = 1;
j = 1;

q_init = 0*pi/180;
theta_init = thets_f(i,j);%45*pi/180;

q    = angle2quat(0,theta_init,phis_f(i,j))';
Do   = quat2dcm(q');
wb   = [0 q_init 0]';
vb   = Do*[v0s(i);0;0];
xj0  = [vb; wb; q];
wj0  = zeros(3,1);
uj0  = [w1s_f(i,j) w2s_f(i,j) d1s_f(i,j) d2s_f(i,j)]';

% Initial Conditions
theta_CI = theta_init;
phi_CI   = phis_f(i,j);
vx_CI    = xj0(1);
vy_CI    = xj0(2);
vz_CI    = xj0(3);

%% Run simulaton Darko
disp('Running the simulation ...')
sim('Darko.slx')

%% Plot 3D trajectory
disp('Plotting results ...')
figure(1)
plot3(pos_x_d, pos_y_d, pos_z_d, 'r');
hold on;
plot3(positions(:,1), positions(:,2), positions(:,3));
grid on;
legend('Setpoint','Measure')
xlabel('North [m]')
ylabel('East [m]')
zlabel('Altitude [m]')