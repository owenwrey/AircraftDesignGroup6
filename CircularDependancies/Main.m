% Circular Dependancies Solver Main equation
% Aircraft Design - Chakraborty
% Group 6
%--------------------------------------------------------------------------
clc; clear; close all

addpath(genpath('Functions')); % lets matlab see all the functions within Functions folder

%% Mission Select
missionToRun = "Strike"; % either 'A2A', or 'Strike', or 'Both' to use the constraining mission
% missionToRun = "A2A";

%% -| Aircraft Struct (iterated variables)|--------------------------------
% All variables that may change with each iteration should be defined here
% (e.g aircraft.cg.x, aircraft.weight.total, aircraft.weight.empty)
aircraft = struct;

% Geometry

% Empty Weight Buildup

% CG and Inertia
aircraft.cg.x = 0;
aircraft.cg.y = 0;
aircraft.cg.z = 0;

% Landing Gear

% Time-Step Mission
aircraft.weight.total = 60e3;
aircraft.weight.fuel = 20e3;
aircraft.weight.empty = 40e3;
aircraft.weight.totalOnLanding = 42e3;

%% -| Aircraft Struct (constant "variables")|------------------------------
% All variables that do not change in the loop but are necessary for
% calculations in the loop (e.g. aircraft.constants.wingLoading, aircraft.weight.tolerance)

% General
aircraft.constants.wingLoading = 112; % [lbf/ft]
aircraft.engine.weight = 3920;
aircraft.engine.thrust = 29000;
aircraft.engine.TSFC = .67;

aircraft.constants.fuelVolume = 3500;
aircraft.fuelSys.VI = 0;                % integral fuel volume, gal
aircraft.fuelSys.VP = aircraft.constants.fuelVolume/2;  % self-sealing tanks volume, gal
aircraft.fuelSys.Nt = 4;

% Geometry
aircraft.fuselage.length   = 48;
aircraft.fuselage.diameter = 5.4;

aircraft.wing.AR          = 8.5;
aircraft.wing.taper_ratio = 0.35;
aircraft.wing.sweep = 30;

aircraft.ht.VolCoeff      = 0.40;
aircraft.ht.AR            = 4.0;
aircraft.ht.TaperRatio    = 0.60;
aircraft.ht.leverArm_frac = 0.45;
aircraft.ht.sweep = 30;

aircraft.vt.VolCoeff      = 0.04;
aircraft.vt.AR            = 1.8;
aircraft.vt.TaperRatio    = 0.70;
aircraft.vt.leverArm_frac = 0.40;
aircraft.vt.twinTail      = true;
aircraft.vt.sweep = 35;

% Empty Weight Buildup
aircraft.constants.limitLoad = 8;
aircraft.constants.ultLoad = 1.5 * aircraft.constants.limitLoad;
aircraft.constants.maxMach = 1.6;

% CG and Inertia
aircraft.engine.cg.x = 45;
aircraft.cockpit.cg.x = 8;

% Landing Gear
aircraft.gear.mg.extendedLength = 40;   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% fixxxxxxxxxxxxxxxxxxxx
aircraft.gear.ng.extendedLength = 40;

% Time-Step Mission
aircraft.weight.tolerance = 250; % GO TO getConfig to uncomment
% -------------------------------------------------------------------------

%% Calculation Loop

exitFlag = false; 
iteration = 0;
iterationMax = 1000;

while( not(exitFlag) && iteration < iterationMax )
    iteration = iteration + 1;

    aircraftOld = aircraft;
    
    %% -| Geometry Updater |-----------------------------------------------
    aircraft = dimensionalize_aircraft(aircraft);
    %----------------------------------------------------------------------


    %%-| Aero Updater |----------------------------------------------------

    %----------------------------------------------------------------------



    %% -| Empty Weight Buildup |-------------------------------------------
    aircraft = EWB(aircraft); 
    %----------------------------------------------------------------------



    %% -| CG and Inertia Calculator |--------------------------------------

    %----------------------------------------------------------------------



    %% -| Landing Gear Updater |-------------------------------------------

    %----------------------------------------------------------------------

    
    %% -| Landing Gear Convergence Check |---------------------------------

    %----------------------------------------------------------------------


    %% -| Fixed MTOW convergence Check |-----------------------------------
    if abs(aircraftOld.weight.total - aircraft.weight.total) > aircraft.weight.tolerance
      continue; % this should go back to the top of the while loop
    end % go on to time-iterated mission model
    %----------------------------------------------------------------------


    %% -| Time Iterated Mission Model |------------------------------------
    aircraft = TIMESTEP_CONVERGENCE_MASTER(aircraft, missionToRun);
    %----------------------------------------------------------------------


    %% -| Converged Solution Check |---------------------------------------
    if abs(aircraftOld.weight.total - aircraft.weight.total) < aircraft.weight.tolerance
      exitFlag = true;
    end
    %----------------------------------------------------------------------

end


%% -| Display Results |----------------------------------------------------
fprintf('\n Converged after %u iterations\n', iteration)
%--------------------------------------------------------------------------

