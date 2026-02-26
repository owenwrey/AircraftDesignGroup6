function aircraft = TIMESTEP_CONVERGENCE_MASTER(aircraft, missionToRun)
% inputs: aircraft (structure)
%         missionToRun: a string that is either 'A2A' or 'Strike' or 'Both'

%% Initialize struct and constants
% Ensure required sub-structs exist
if ~isfield(aircraft, 'constants') || ~isstruct(aircraft.constants)
    aircraft.constants = struct();
end

if ~isfield(aircraft, 'weight') || ~isstruct(aircraft.weight)
    aircraft.weight = struct();
end

% Only set defaults if they do not already exist
if ~isfield(aircraft.constants, 'wingLoading')
    aircraft.constants.wingLoading = 115;
    warning("Used hardcoded wingLoading")
end

if ~isfield(aircraft.weight, 'tolerance')
    aircraft.weight.tolerance = 0.06;
    warning("Used hardcoded weight.tolerance")
end

%% Run Mission

missionToRun = string(missionToRun);

% switch looks at which mission is set to run, runs it,
% and replaces aircraft with its aircraft output
switch missionToRun

    case "A2A"
        aircraft = A2A(aircraft);

    case "Strike"
        aircraft = Strike(aircraft);

    case "Both"
        Results = struct;
        Results.A2A = A2A(aircraft);       
        Results.Strike = Strike(aircraft);

        if Results.A2A.weight.total > Results.Strike.weight.total
            aircraft.constrainingMission = "A2A";
            aircraft = Results.A2A;
        elseif Results.Strike.weight.total > Results.A2A.weight.total
            aircraft.constrainingMission = "Strike";
            aircraft = Results.Strike;
        end

    otherwise
        error("Unknown mission type: %s", missionToRun)

end % switch


end % function