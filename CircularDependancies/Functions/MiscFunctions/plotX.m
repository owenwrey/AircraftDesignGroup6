function plotX(aircraft)
% ylines x coords of components

field = fieldnames(aircraft)';
blacklist = ["flightcond","cg","constants","aero","enginesystems","weight","fuelSys","gear","TimeStepTable"];

figure
hold on

for i = 1:numel(field)
    if ~ismember(blacklist,field{i})
        scatter(aircraft.(field{i}).cg.x,5,DisplayName=field{i})
    end
end

xline(aircraft.cg.x,DisplayName="CG")

scatter(aircraft.gear.mg.cg.x,5,DisplayName="mg")
scatter(aircraft.gear.ng.cg.x,5,DisplayName="ng")
legend
plot([1 48],[0 0])
xline(48)
end