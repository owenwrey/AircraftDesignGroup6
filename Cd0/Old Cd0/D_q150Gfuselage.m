function [D_q] = D_q150Gfuselage(M)
%D_Q150Gfuselage Calculates Drag over dynamic pressure for 150 gallon wing tank
%   Inputs: M - Mach number
%   Outputs: D/q
%   taken from P.426 in Raymer textbook.
arguments (Input)
    M double
end

arguments (Output)
    D_q
end
Data = [0.425,0.297
0.45,0.2975
0.475,0.2995
0.5,0.302
0.525,0.3038
0.55,0.3044
0.575,0.3041
0.6,0.3035
0.625,0.3029
0.65,0.3027
0.675,0.3031
0.7,0.3041
0.725,0.3045
0.75,0.3022
0.775,0.2952
0.8,0.2928
0.825,0.3185
0.85,0.3751
0.875,0.5354
0.9,0.7512
0.925,1.075
0.95,1.4231
0.975,1.6084
1,1.69];

func =griddedInterpolant(Data(:,1),Data(:,2));

D_q = func(M);


end