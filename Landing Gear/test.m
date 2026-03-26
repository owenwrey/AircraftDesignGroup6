% landing gear
clc; clear; close all;
aircraft.cg.x = 25; % ft
aircraft.cg.z = 8; % ft
aircraft.cg.y = 0; % ft
aircraft.wing.span = 40; % ft
aircraft.fuselage.radius = 4; % ft 
aircraft.fuselage.length = 50; %ft

Aircraft = landingGear(aircraft);