function [D_q] = D_q150Gwing(M)
%D_Q150Gwing Calculates Drag over dynamic pressure for 150 gallon wing tank
%   Inputs: M - Mach number
%   Outputs: D/q
%   taken from P.426 in Raymer textbook.
arguments (Input)
    M double
end

arguments (Output)
    D_q
end
Data = [0.425,0.3339
0.45,0.334
0.475,0.3346
0.5,0.3361
0.525,0.3384
0.55,0.3407
0.575,0.3419
0.6,0.342
0.625,0.3415
0.65,0.3409
0.675,0.3408
0.7,0.3422
0.725,0.3434
0.75,0.3407
0.775,0.3307
0.8,0.325
0.825,0.3531
0.85,0.4461
0.875,0.623
0.9,0.8857
0.925,1.2686
0.95,1.5664
0.975,1.7395];

func =griddedInterpolant(Data(:,1),Data(:,2));

D_q = func(M);


end