% Airfoil area calculator
%
%
%--------------------------------------------------------------------------
clc; clear; close all

data = importdata("naca64206.txt");

airfoil = polyshape(data.data); % chord length of 1

Area = area(airfoil)

plot(airfoil)
axis equal
