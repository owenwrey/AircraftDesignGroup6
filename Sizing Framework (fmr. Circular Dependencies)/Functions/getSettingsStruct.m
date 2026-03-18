function sett = getSettingsStruct()
% Settings for plots and which other programs to run after convergence

% FORMAT:
% sett.SCRIPT_NAME.run ---> determines if SCRIPT_NAME will be run
% sett.SCRIPT_NAME.plot ---> determines if SCRIPT_NAME plots will be displayed
% sett.SCRIPT_NAME.print ---> determines if SCRIPT_NAME will print to command window
sett = struct();

% Add any scripts I missed!
%% CG Envelope
sett.CGenvelope.run = 1;
sett.CGenvelope.plot = 0;
sett.CGenvelope.print = 0;

%% Time-iterated Mission Model
sett.TimeStep.plot.fuelBurn = 0;
sett.TimeStep.plot.weightFraction = 0;
sett.TimeStep.print = 0;

%% Excess Power
% Time to Climb settings
sett.TimeToClimb.run = 1;
sett.TimeToClimb.plot = 0;
sett.TimeToClimb.print = 0;

% Alt-Speed Envelope settings
sett.AltSpeedEnvelope.run = 1;
sett.AltSpeedEnvelope.plot = 0;
sett.AltSpeedEnvelope.print = 0;

% Turn Performance settings
sett.TurnPerformace.run = 1;
sett.TurnPerformace.plot = 0;
sett.TurnPerformace.print = 0;

%% V-n Diagram
sett.VnDiagram.run = 1;
sett.VnDiagram.plot = 0;
sett.VnDiagram.print = 0;

end