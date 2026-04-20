function CGenvelope(aircraft, Title)
%CGENVELOPE calculates and plots cg envelope for ordinance drop mission

figure

drop = true;
time = aircraft.TimeStepTable.Time;
dfuel = (aircraft.TimeStepTable.dFuel);
weight = aircraft.TimeStepTable.Weight;

 % bodge fix


blacklist = ["fuelLength","totalfuelvolume","flightcond","cg","constants","aero","enginesystems","weight","fuelSys","gear","TimeStepTable","inertia"];

CGxloc = zeros(1,length(time));

sumVolume = aircraft.fuselagetank1.volume + aircraft.fuselagetank2.volume + aircraft.wingtanks.volume;

FTfact = aircraft.fuselagetank1.volume./sumVolume; % fuel tank factor
WTfact = aircraft.wingtanks.volume./sumVolume; % wing tank factor

for j = 1 : length(time)

    aircraft.fuselagetank1.weight = aircraft.fuselagetank1.weight - dfuel(j)*FTfact;
    aircraft.fuselagetank2.weight = aircraft.fuselagetank2.weight - dfuel(j)*FTfact;
    aircraft.wingtanks.weight = aircraft.wingtanks.weight - dfuel(j)*WTfact;

   

    if drop == true && j == 42

        aircraft.ordinance.weight = 0;

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

figure(10)

forwardLimit = aircraft.cg.x-aircraft.wing.MAC*(0.5);

statMarg = 0.15;
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