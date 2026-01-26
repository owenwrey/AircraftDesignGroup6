%% main_cd0.m
% Mission CD0 calculator with subsonic/supersonic switch

clear; close all; clc;

%% Params (constants)
params.mu   = 1.789e-5;     % dynamic viscosity [kg/(m*s)]
params.Sref = 578.846;      % reference area [ft^2]
params.k_default = 0.50e-5; % roughness [ft] polished sheet metal

%% Mission inputs
alt_ft = [ ...
    0
    0
    25000
    25000
    25000
    0
    25000
    25000
    25000
    0
    0 ];

M = [0 0 0.8 1.1 0.8 0.85 0.8 0.94 0.8 0.42 0];

%% Build components
comp = buildComponents(params.k_default,M);

%% Compute CD0 vector
Cd0_total = computeCd0Mission(alt_ft, M, comp, params);

disp("CD0 total per segment:")
disp(Cd0_total)
