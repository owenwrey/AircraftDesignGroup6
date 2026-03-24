% Circular Dependancies Solver Main equation
% Aircraft Design - Chakraborty
% Group 6
%--------------------------------------------------------------------------
clc; clear; close all

total = tic;

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
aircraft.cg.x = 20;
aircraft.cg.y = 0;
aircraft.cg.z = 10;

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
aircraft.fuselage.diameter = 6;
aircraft.fuselage.cg.x = 24;
aircraft.fuselage.cg.y = 0;
aircraft.fuselage.cg.z = 10;

aircraft.wing.AR          = 4;
aircraft.wing.taper_ratio = 0.25;
aircraft.wing.sweep = 30;
aircraft.wing.T2C = .055;   % thickness to chord
aircraft.wing.l = 50;
aircraft.wing.x_c = .24;

aircraft.ht.VolCoeff      = 0.40;
aircraft.ht.AR            = 4.0;
aircraft.ht.TaperRatio    = 0.40;
aircraft.ht.leverArm_frac = 0.35;
aircraft.ht.sweep = 30;
aircraft.ht.T2C = .05;
aircraft.ht.x_c = .24;

aircraft.vt.VolCoeff      = 0.04;
aircraft.vt.AR            = 1.8;
aircraft.vt.TaperRatio    = 0.30;
aircraft.vt.leverArm_frac = 0.30;
aircraft.vt.twinTail      = true;
aircraft.vt.sweep = 35;
aircraft.vt.T2C = .05;
aircraft.vt.x_c = .24;

aircraft.ordinance.weight = 4380;

aircraft.avionics.weight = 2500;

% Empty Weight Buildup
aircraft.constants.limitLoad = 8;
aircraft.constants.ultLoad = 1.5 * aircraft.constants.limitLoad;
aircraft.constants.maxMach = 1.6;

% CG and Inertia
aircraft.engine.cg.x = 45;
aircraft.cockpit.cg.x = 8;
aircraft.cockpit.weight = 300; % undefined? %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Landing Gear
aircraft.gear.mg.extendedLength = 40;   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% fixxxxxxxxxxxxxxxxxxxx
aircraft.gear.ng.extendedLength = 40;
aircraft.gear.mg.x = 40;
aircraft.gear.ng.x = 8;
aircraft.gear.mg.height = 5;
aircraft.gear.ng.height = 5; 

% Time-Step Mission
aircraft.weight.tolerance = 5; % GO TO getConfig to uncomment
aircraft.cg.tolerance = 1/12;
% -------------------------------------------------------------------------

%% Calculation Loop

exitFlag = false; 
iteration = 0;
iterationMax = 100;



while( not(exitFlag) && iteration <= iterationMax )
    iteration = iteration + 1;
    fprintf("   Iteration: %u\n", iteration);
    aircraftOld = aircraft;


    tic
    fprintf("Geometry\n")
    %% -| Geometry Updater |-----------------------------------------------
    aircraft = dimensionalize_aircraft(aircraft);
    %----------------------------------------------------------------------
    toc

        tic
    fprintf("aero\n")
    %% -| Aero Updater |----------------------------------------------------
    aircraft = aeroupdater(aircraft);
    %----------------------------------------------------------------------
    toc

        tic
    fprintf("EWB\n")
    %% -| Empty Weight Buildup |-------------------------------------------
    aircraft = EWB(aircraft); 
    %----------------------------------------------------------------------
    toc

        tic
    fprintf("CGInertia\n")
    %% -| CG and Inertia Calculator |--------------------------------------
     aircraft = CgInertiaCalc(aircraft);
    %----------------------------------------------------------------------
    toc

        tic
    fprintf("Landing Gear\n")
    %% -| Landing Gear Updater |-------------------------------------------
    aircraft = landingGear(aircraft);
    %----------------------------------------------------------------------
    toc
    fprintf(" \n\n")

    %% -| Landing Gear Convergence Check |---------------------------------

    %----------------------------------------------------------------------

  

    %% -| Fixed MTOW Convergence Check |-----------------------------------
    if abs(aircraftOld.weight.total - aircraft.weight.total) > aircraft.weight.tolerance

      continue; % this should go back to the top of the while loop
     
    end % go on to time-iterated mission model
    %----------------------------------------------------------------------

    %% -| Time Iterated Mission Model |------------------------------------
    aircraft = TIMESTEP_CONVERGENCE_MASTER(aircraft, missionToRun);
    %----------------------------------------------------------------------

    %% -| Converged Solution Check |---------------------------------------
    if abs(aircraftOld.weight.total - aircraft.weight.total) < aircraft.weight.tolerance
        if abs(aircraftOld.cg.x - aircraft.cg.x) < aircraft.cg.tolerance
      exitFlag = true;
        end
    end
    %----------------------------------------------------------------------
    
    aircraft.constants.fuelVolume = aircraft.weight.fuel/6.7;
    aircraft.fuelSys.VP = aircraft.constants.fuelVolume/2;  % self-sealing tanks volume, gal

end


% ignore this, just helps with report writing
% cell2mat(struct2cell(aircraft.engine.structures))
% cell2mat(struct2cell(aircraft.weight.misc))
% enginePerc = cell2mat(struct2cell(aircraft.engine.structures))./aircraft.weight.total
% miscPerc = cell2mat(struct2cell(aircraft.weight.misc))./aircraft.weight.total

%% -| Display Results |----------------------------------------------------
fprintf("\n Converged after %u iterations\n", iteration)
fprintf("  MTOW: %.0f lb\n", aircraft.weight.total)
fprintf("  EOW: %.0f lb\n", aircraft.weight.empty)
fprintf("  CG_x: %.3f ft\n", aircraft.cg.x)
fprintf("  CG_y: %.3f ft\n", aircraft.cg.y)
fprintf("  CG_z: %.3f ft\n", aircraft.cg.z)
%--------------------------------------------------------------------------

% plot cg envelope

if false
CGenvelope(aircraft, "Strike, No Drop")
end

runtime = toc(total);

fprintf("\n\n The total runtime is : " + runtime)