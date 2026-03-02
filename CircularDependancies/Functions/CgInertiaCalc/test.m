% Test script for CG function
% Aircraft design group 6
% Chakraborty - 2026
clc; clear; close all

aircraft.wing.x = 10;
aircraft.wing.y = 5;
aircraft.wing.z = 3;
aircraft.wing.weight = 100;

aircraft.fuel.x = 18;
aircraft.fuel.y = 9;
aircraft.fuel.z = 3;
aircraft.fuel.weight = 300;

aircraft.constants.cd0 = 0.05;


out = CgInertiaCalc(aircraft)
