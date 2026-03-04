function CGenvelope(aircraft, Title)
%CGENVELOPE calculates and plots cg envelope for ordinance drop mission

figure

time = aircraft.TimeStepTable.Time;
dfuel = aircraft.TimeStepTable.dFuel;
weight = aircraft.TimeStepTable.Weight;

halfFuel = aircraft.weight.fuel/2;


blacklist = ["flightcond","cg","constants","aero","enginesystems","weight","fuelSys","gear","TimeStepTable"];

CGxloc = zeros(1,length(time));

for j = 1 : length(time)

    aircraft.fuel2.weight = aircraft.fuel2.weight - dfuel(j);

    if aircraft.fuel2.weight <= 0
        aircraft.fuel1.weight = aircraft.fuel1.weight + aircraft.fuel2.weight;
        aircraft.fuel2.weight = 0;
    end

    % initialize variables
xsum = 0;
ysum = 0;
zsum = 0;
weightsum = 0;

field = fieldnames(aircraft)';

for i = 1:numel(field)

    
    if ~ismember(blacklist,field{i})

        
        % sum weights times locations
        xsum = aircraft.(field{i}).weight.*aircraft.(field{i}).cg.x + xsum;
        
        % sum weight
        weightsum = aircraft.(field{i}).weight + weightsum;
       
    end
end

%landing gear

    % NG sum weights times locations
        xsum = aircraft.gear.mg.weight.*aircraft.gear.mg.cg.x + xsum;
        
        % sum weight
        weightsum = aircraft.gear.mg.weight + weightsum;

   % MG sum weights times locations
        xsum = aircraft.gear.ng.weight.*aircraft.gear.ng.cg.x + xsum;

        % sum weight
        weightsum = aircraft.gear.ng.weight + weightsum;
        

% cg  & add to output
CGxloc(j) = xsum./weightsum;



end

forwardLimit = aircraft.cg.x-aircraft.wing.MAC/2;

statMarg = 0.05;
aftLimit = aircraft.cg.x + aircraft.wing.MAC*statMarg;
plot(CGxloc,weight)
hold on
xline(forwardLimit)
xline(aftLimit)

title("Mission: " + Title)
xlabel("CG fuselage station (ft.)")
ylabel('Gross Weight (lbs)');
legend('CG Location', 'Forward Limit', 'Aft Limit');
grid on;

end