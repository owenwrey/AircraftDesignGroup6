function aircraft = TIMESTEP_CONVERGENCE_MASTER(aircraft, missionToRun)
% inputs: aircraft (structure)
%         missionToRun: a string that is either 'A2A' or 'Strike' or 'Both'

Results = struct;

missionToRun = string(missionToRun);

% switch looks at which mission is set to run, runs it,'
% and replaces aircraft with its aircraft output
switch missionToRun

    case "A2A"
        Results.A2A = A2A(aircraft);
        aircraft = Results.A2A;

    case "Strike"
        Results.Strike = Strike(aircraft);
        aircraft = Results.Strike;

    case "Both"
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