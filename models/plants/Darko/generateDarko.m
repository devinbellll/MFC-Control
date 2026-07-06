function drone = generateDarko( )
%GENERATEMAVION generate Dark-0 drone data structure
%
%  this function generates mavion physical parameters structure for use 
%  with the Phi-theory toolbox.

%% initialize Dark-0 drone structure

                          % 1 => MAVion
drone_specifications = 2; % 2 => Darko

if(drone_specifications == 1)
% MAVion
drone = struct('MASS',               0.45,      ... % ok: 0.438
               'INERTIA',            eye(3),     ... % below
               'G',                  9.81,       ... % ok
               'P_P1_CG',            [1; 1; 1],  ... % below
               'P_P2_CG',            [1; 1; 1],  ... % below
               'P_A1_CG',            [1; 1; 1],  ... % below
               'P_A2_CG',            [1; 1; 1],  ... % below
               'INERTIA_PROP_X',     3.46e-6,    ... % ok
               'INERTIA_PROP_N',     0,          ... % ok
               'PHI',                eye(6),     ... % below
               'RHO',                1.225,      ... % ok
               'WET_SURFACE',        0.0441,     ... % ok
               'DRY_SURFACE',        0,          ... % ok
               'PHI_n',              0,          ... % ok
               'CHORD',              0.21,       ... % ok
               'WINGSPAN',           0.42,       ... % ok
               'PROP_RADIUS',        0.105,      ... % ok
               'ELEVON_MEFFICIENCY', [0;0.66;0], ... % ok 0.8
               'ELEVON_FEFFICIENCY', [0;0.33;0], ... % ok 0.43
               'PROP_KP',            4.48e-06,...    % ok: 1.75*4.48e-06
               'THICKNESS',          0.02    ,   ... % ok
               'PROP_KM',            2.400e-7     ); % ok
           %% acording to back-of-the-envelope computations
drone.INERTIA = diag([0.0036 0.0036 0.0072]);

end
if(drone_specifications == 2)
%Darko
drone = struct('MASS',               0.492,      ... % ok: 0.438
               'INERTIA',            eye(3),     ... % below
               'G',                  9.81,       ... % ok
               'P_P1_CG',            [1; 1; 1],  ... % below
               'P_P2_CG',            [1; 1; 1],  ... % below
               'P_A1_CG',            [1; 1; 1],  ... % below
               'P_A2_CG',            [1; 1; 1],  ... % below
               'INERTIA_PROP_X',     5.1116e-6,    ... % ok
               'INERTIA_PROP_N',     0,          ... % ok
               'PHI',                eye(6),     ... % below
               'RHO',                1.225,      ... % ok
               'WET_SURFACE',        0.0743,     ... % ok
               'DRY_SURFACE',        0,          ... % ok
               'PHI_n',              0,          ... % ok
               'CHORD',              0.13,       ... % ok
               'WINGSPAN',           0.55,       ... % ok
               'PROP_RADIUS',        0.125,      ... % ok
               'ELEVON_MEFFICIENCY', [0;0.93;0], ... % ok
               'ELEVON_FEFFICIENCY', [0;0.48;0], ... % ok
               'PROP_KP',            5.13e-6,...    % ok: 1.75*4.48e-06
               'THICKNESS',          0.02   ,   ... % ok
               'PROP_KM',            2.640e-7     ); % ok
           %% acording to inertia identificayiton (python + IMU measurements)
           drone.INERTIA = diag([0.007018 0.002785 0.00606]);
end

%% according to thin airfoil phi-theory (see [1])
Cd0 = 0.025;
Cy0 = 0.1; 
% Cd0 = 0.1; %MAVION 
% Cy0 = 0.1; %MAVION 

dR = -0.1*drone.CHORD;
% PHI_fv = diag([Cd0; Cy0; (2*pi+Cd0)]);
% PHI_mv = [0 0 0 ; 0 0 -1/drone.CHORD*dR*(2*pi+Cd0); 0 1/drone.WINGSPAN*dR*Cy0 0];
% %PHI_mw = 1/2*diag([0.5 0.5 0.5]); %Cl Cm Cn MAVION
% %PHI_mw = 1/2*diag([0.47 0.54 0.52]); %Cl Cm Cn
% PHI_mw = 1/2*([0.2758 0 0.1194; 0 0.8825 0; 0.0665 0 0.0044]); %Cl Cm Cn
% drone.PHI = [ PHI_fv PHI_mv; PHI_mv PHI_mw ];

% Revisited Phi-theory with Prandtl's lifting-line theory
AR    = (drone.WINGSPAN^2)/drone.WET_SURFACE;
%Pr_c  = 2*pi/(1+2/(0.95*AR));
die_c = pi*AR/(1+sqrt(1+(AR^2/4)));
PHI_fv = diag([Cd0; Cy0; (die_c+Cd0)]);
PHI_mv = [0 0 0 ; 0 0 -1/drone.CHORD*dR*(die_c+Cd0); 0 1/drone.WINGSPAN*dR*Cy0 0];
%PHI_mw = 1/2*diag([0.5 0.5 0.5]); %Cl Cm Cn MAVION
%PHI_mw = 1/2*diag([0.47 0.54 0.52]); %Cl Cm Cn
PHI_mw = 1/2*([0.2792 0 0.1145; 0 1.2715 0; 0.081 0 0.0039]); %Cl Cm Cn 10% static margin
drone.PHI = [PHI_fv PHI_mv; PHI_mv PHI_mw];

%% geometric parameters
% drone.P_P1_CG = [0.15; -0.5; 0]*drone.CHORD; % 10% static margin MAVion
% drone.P_P2_CG = [0.15;  0.5; 0]*drone.CHORD; % 10% static margin MAVion
% drone.P_A1_CG = [ 0.0; -0.5; 0]*drone.CHORD; %                   MAVion 
% drone.P_A2_CG = [ 0.0;  0.5; 0]*drone.CHORD; %                   MAVion
drone.P_P1_CG = [0.065; -0.155; 0];             % 10% static margin Dark0
drone.P_P2_CG = [0.065;  0.155; 0];             % 10% static margin Dark0
drone.P_A1_CG = [ 0.0; -0.155; 0];             %                   Dark0
drone.P_A2_CG = [ 0.0;  0.155; 0];             %                   Dark0

end