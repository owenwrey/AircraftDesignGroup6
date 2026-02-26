% Circular Dependancies Solver Main equation
% Aircraft Design - Chakraborty
% Group 6
%--------------------------------------------------------------------------
clc; clear; close all

addpath(genpath('Functions')); % lets matlab see all the functions within Functions folder

%% Mission Select
missionToRun = "Strike"; % either 'A2A', or 'Strike', or 'Both' to use the constraining mission
% missionToRun = "A2A";

%% -| Aircraft Struct |----------------------------------------------------
% All variable information related to the aircraft should be
% stored/accesible in this struct.

aircraft.cg.x = 0;
aircraft.cg.y = 0;
aircraft.cg.z = 0;

% This just makes sure we are getting the tolerance from the same function
% (getConfig) the time-step is. If you need to change it, change it there.
tempConfig = getConfig();
aircraft.weight.tolerance = tempConfig.weightTolerance;
clearvars tempConfig;

%example component



% -------------------------------------------------------------------------


%% Calculation Loop

exitFlag = false; 
iteration = 0;
iterationMax = 1000;

while( not(exitFlag) && iteration < iterationMax )
    iteration = iteration + 1;

    aircraftOld = aircraft;
    
    %-| Geometry Updater |-------------------------------------------------
    aircraft = dimensionalize_aircraft(aircraft);
    %----------------------------------------------------------------------


    %-| Aero Updater |-----------------------------------------------------

    %----------------------------------------------------------------------



    %-| Empty Weight Buildup |---------------------------------------------
    aircraft = EWB(aircraft); % idk if these are the right i/o's - Owen
    %----------------------------------------------------------------------



    %-| CG and Inertia Calculator |----------------------------------------

    %----------------------------------------------------------------------



    %-| Landing Gear Updater |---------------------------------------------

    %----------------------------------------------------------------------

    
    %-| Landing Gear Convergence Check |-----------------------------------

    %----------------------------------------------------------------------


    %-| Fixed MTOW convergence Check |-------------------------------------
    if abs(aircraftOld.weight.total - aircraft.weight.total) > aircraft.weight.tolerance
      continue; % this should go back to the top of the while loop
    end % go on to time-iterated mission model
    %----------------------------------------------------------------------


    %-| Time Iterated Mission Model |--------------------------------------
    aircraft = TIMESTEP_CONVERGENCE_MASTER(aircraft, missionToRun);
    %----------------------------------------------------------------------


    %-| Converged Solution Check |-----------------------------------------
    if abs(aircraftOld.weight.total - aircraft.weight.total) < aircraft.weight.tolerance
      exitFlag = true;
    end
    %----------------------------------------------------------------------

end


%-| Display Results |------------------------------------------------------
fprintf('Converged after %u iterations\n', iteration)
%--------------------------------------------------------------------------
