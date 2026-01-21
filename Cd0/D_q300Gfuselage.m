function [D_q] = D_q300Gfuselage(M)
%D_Q300Gfuselage Calculates Drag over dynamic pressure for 300 gallon fuselage tank
%   Inputs: M - Mach number
%   Outputs: D/q
%   taken from P.426 in Raymer textbook.
arguments (Input)
    M double
end

arguments (Output)
    D_q
end
Data = [0.425,0.4169
0.45,0.4121
0.475,0.4098
0.5,0.4092
0.525,0.4095
0.55,0.4099
0.575,0.4096
0.6,0.4092
0.625,0.4093
0.65,0.4107
0.675,0.4125
0.7,0.4123
0.725,0.4075
0.75,0.3997
0.775,0.3963
0.8,0.4049
0.825,0.4376
0.85,0.5272
0.875,0.7131
0.9,1.0288
0.925,1.4542
0.95,1.8507
0.975,2.0863
1,2.1873
];

func = griddedInterpolant(Data(:,1),Data(:,2));

D_q = func(M);

end