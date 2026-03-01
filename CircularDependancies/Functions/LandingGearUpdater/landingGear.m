function [Aircraft] = landingGear(Aircraft)

tipBackAngle = 20;

WB = 0.8.* Aircraft.fuselage.length;

J = WB./10; % main gear moment arm
K = WB.*(9/10); % front gear moment arm

Aircraft.mainGear.x = Aircraft.cg.x + J;
Aircraft.frontGear.x = Aircraft.cg.x - K;

Aircraft.mainGear.height = (Aircraft.fuselage.length - Aircraft.mainGear.x).*tand(tipBackAngle);
Aircraft.frontGear.height = Aircraft.mainGear.height; 

h_wingtip = Aircraft.cg.z;  
b_2 = Aircraft.wing.span/2; 
Aircraft.mainGear.width = 2*(b_2 - (h_wingtip-0.5))./(tand(5)); % 5 degree roll


end