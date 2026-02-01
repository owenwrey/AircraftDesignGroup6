function [CD0] = CD0F18(M)
%D_Q2kBfuselage Calculates CD0 for the f18 at a given mach number
%   Inputs: M - Mach number
%   Outputs: D/q
%   taken from P.426 in Raymer textbook.
arguments (Input)
    M double
end

arguments (Output)
    CD0
end
Data = [0.6,0.0251
0.7,0.0249
0.8,0.0253
0.9,0.0283
1,0.0437
1.1,0.057
1.2,0.0586
1.3,0.0585
1.4,0.0583
1.5,0.0581
1.6,0.058
];

func =griddedInterpolant(Data(:,1),Data(:,2));

CD0 = func(M);


end