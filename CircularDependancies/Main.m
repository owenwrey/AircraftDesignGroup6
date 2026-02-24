% Circular Dependancies Solver Main equation
% Aircraft Design - Chakraborty
% Group 6
%--------------------------------------------------------------------------
clc; clear; close all

addpath(genpath('Functions')); % lets matlab see all the functions within Functions folder

%% Mission Select
missionToRun = 'Strike'; % either 'A2A', or 'Strike', or 'Both' to use the constraining mission

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

    %----------------------------------------------------------------------



    %-| Aero Updater |-----------------------------------------------------

    %----------------------------------------------------------------------



    %-| Empty Weight Buildup |---------------------------------------------

    %----------------------------------------------------------------------



    %-| CG and Inertia Calculator |----------------------------------------

    %----------------------------------------------------------------------



    %-| Landing Gear Updater |---------------------------------------------

    %----------------------------------------------------------------------


    %-| Fixed MTOW convergence Check |-------------------------------------
    % type "help continue" to see how to send while loop back to top

    if abs(aircraftOld.weight.total - aircraft.weight.total) > aircraft.weight.tolerance
      continue; % this should go back to the top of the while loop
    end % go on to time-iterated mission model


    %----------------------------------------------------------------------



    %-| Time Iterated Mission Model |--------------------------------------
    % run time step here
    aircraft = TIMESTEP_CONVERGENCE_MASTER(aircraft, missionToRun);
    %----------------------------------------------------------------------


    %-| Converged Solution Check |-----------------------------------------
    if abs(aircraftOld.weight.total - aircraft.weight.total) < aircraft.weight.tolerance
      exitFlag = true;
    end
    %----------------------------------------------------------------------

end


%-| Display Results |------------------------------------------------------
fprintf('Converged!!!\n')
%--------------------------------------------------------------------------