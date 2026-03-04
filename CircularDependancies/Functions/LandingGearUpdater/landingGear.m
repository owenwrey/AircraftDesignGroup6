function [Aircraft] = landingGear(Aircraft)

tipBackAngle = 20;

WB = 0.8.* Aircraft.fuselage.length;

J = WB./10; % main gear moment arm
K = WB.*(9/10); % front gear moment arm

Aircraft.gear.mg.x = Aircraft.cg.x + J;
Aircraft.gear.ng.x = Aircraft.cg.x - K;

Aircraft.gear.mg.height = (Aircraft.fuselage.length - Aircraft.gear.mg.x).*tand(tipBackAngle);
Aircraft.gear.ng.height = Aircraft.gear.mg.height; 

h_wingtip = Aircraft.cg.z;  
b_2 = Aircraft.wing.span/2; 
Aircraft.gear.mg.width = 2*(b_2 - (h_wingtip-0.5))./(tand(5)); % 5 degree roll


end