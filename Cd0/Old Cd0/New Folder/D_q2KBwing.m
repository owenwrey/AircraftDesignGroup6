function [D_q] = D_q2KBwing(M)
%D_Q2kBwing Calculates Drag over dynamic pressure for 2000 lb bomb on wing
%   Inputs: M - Mach number
%   Outputs: D/q
%   taken from P.426 in Raymer textbook.
arguments (Input)
    M double
end

arguments (Output)
    D_q
end
Data = [0.4,0.1455
0.425,0.1452
0.45,0.146
0.475,0.1477
0.5,0.1501
0.525,0.1531
0.55,0.1562
0.575,0.1595
0.6,0.1625
0.625,0.1651
0.65,0.1671
0.675,0.1682
0.7,0.1683
0.725,0.167
0.75,0.1647
0.775,0.1642
0.8,0.1695
0.825,0.185
0.85,0.22
0.875,0.2894
0.9,0.4117
0.925,0.5835
0.95,0.7766
];

func =griddedInterpolant(Data(:,1),Data(:,2));

D_q = func(M);


end