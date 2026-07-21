clc
clear all
close all

addpath('functions');   % controller math: mfc_siso + the estimator kernels
addpath('blocks');      % matlab.System blocks: mfc_siso_core + the stage blocks
addpath('library');     % Simulink library and its builder
addpath('examples');    % val_mfc, val_mfc_composed
addpath('tests');
addpath('models');
addpath('scripts');
addpath('models/plants/Darko/');

%%
X=1;
Y=2;
Z=3;

load_parameters;

