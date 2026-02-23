function [airOut] = CgInertiaCalc(airIn)
%CGINERTIACALC Calculates CG and products of inertia of Aircraft Struct.
%   Detailed explanation will go here

%% main loop

xsum = 0;
ysum = 0;
zsum = 0;

weightsum = 0;

for field = fieldnames(airIn)'

xsum = airIn.field


end

%% outputs


end