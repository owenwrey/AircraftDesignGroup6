function [D_q] = D_q2KBfuselage(M)
%D_Q2kBfuselage Calculates Drag over dynamic pressure for 2000 lb bomb on
%fuselage
%   Inputs: M - Mach number
%   Outputs: D/q
%   taken from P.426 in Raymer textbook.
arguments (Input)
    M double
end

arguments (Output)
    D_q
end
Data = [0.4,0.1537
0.425,0.155
0.45,0.1559
0.475,0.1567
0.5,0.1573
0.525,0.1579
0.55,0.1584
0.575,0.1591
0.6,0.16
0.625,0.1611
0.65,0.1625
0.675,0.1644
0.7,0.1666
0.725,0.1688
0.75,0.1706
0.775,0.1714
0.8,0.1713
0.825,0.186
0.85,0.2477
0.875,0.3739
0.9,0.5611
0.925,0.7785
0.95,0.9835
];

func =griddedInterpolant(Data(:,1),Data(:,2));

D_q = func(M);


end