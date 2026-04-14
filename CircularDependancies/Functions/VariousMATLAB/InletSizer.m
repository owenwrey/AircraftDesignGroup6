% Inlet Sizer
% supersonic inlet, 2 shock design. pg 290, fig 10.8
% inlet will be variable geometry, have a suck-in door

% Engine front face area: 1075.62 in^2
% Design Mach: .8
% Design Mass flow: 170 lb/s
% Takeoff mass flow: 323 lb/s

clear;clc;

A_engine = 1075.62;     % engine face area, in^2

M_inf = .8;             % design mach
M_engine = .5;          % mach at engine face
M_inlet = (M_inf-M_engine)/2 + M_engine;

ratio_inf = 1.038;
ratio_engine = 1/M_engine*((1+.2*M_engine^2)/1.2)^3;
ratio_inlet = 1/M_inlet*((1+.2*M_inlet^2)/1.2)^3;

inlet_engine = ratio_inlet/ratio_engine;

A_inlet = inlet_engine*A_engine;
A_inlet = A_inlet/144;

size = sqrt(A_inlet);

fprintf("   The inlet area is %.2f ft^2 \n", A_inlet)
