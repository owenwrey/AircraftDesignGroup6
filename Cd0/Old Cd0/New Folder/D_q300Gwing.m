function [D_q] = D_q300Gwing(M)
%D_Q300Gwing Calculates Drag over dynamic pressure for 300 gallon wing tank
%   Inputs: M - Mach number
%   Outputs: D/q
%   taken from P.426 in Raymer textbook.
arguments (Input)
    M double
end

arguments (Output)
    D_q
end
Data = [0.425,0.4858 
0.45,0.4882
0.475,0.4892
0.5,0.4893
0.525,0.4889
0.55,0.4884
0.575,0.4884
0.6,0.4893
0.625,0.4911
0.65,0.4931
0.675,0.4939
0.7,0.4922
0.725,0.4881
0.75,0.4841
0.775,0.4834
0.8,0.4888
0.825,0.5175
0.85,0.617
0.875,0.8176
0.9,1.2042
0.925,1.7948
0.95,2.2368
0.975,2.4751];

func =griddedInterpolant(Data(:,1),Data(:,2));

D_q = func(M);


end