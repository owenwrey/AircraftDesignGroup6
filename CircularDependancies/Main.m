% Circular Dependancies Solver Main equation
% Aircraft Design - Chakraborty
% Group 6
%--------------------------------------------------------------------------
clc; clear; close all

addpath(genpath('Functions')); % lets matlab see all the functions within Functions folder


%% -| Aircraft Struct |----------------------------------------------------
% All variable information related to the aircraft should be
% stored/accesible in this struct.
aircraft = struct;

components = {Wing, Htail, Vtail, fuselage, Gear, Crew, } 
variables




%example component



% -------------------------------------------------------------------------


%% calculation loop

exitFlag = false; 
iteration = 0;
iterationMax = 1000;

while( not(exitFlag) && iteration < iterationMax )
    iteration = iteration + 1;
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

    % Something like this:
    % if abs(aircraftOld.weight.total - aircraft.weight.total) < weightTol
    %   go to next step of loop (time iterated mission model)
    % else
    %   go back to top of loop
    % end


    %----------------------------------------------------------------------



    %-| Time Iterated Mission Model |--------------------------------------
    % run time step here
    %----------------------------------------------------------------------

    

    %-| Converged Solution Check |-----------------------------------------

    %----------------------------------------------------------------------

end


%-| Display Results |------------------------------------------------------

%--------------------------------------------------------------------------