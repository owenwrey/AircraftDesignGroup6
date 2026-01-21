function [D_q] = D_qClusterRack(M)
%D_QAim9 Calculates Drag over dynamic pressure for A multiple bomb cluster
%rack.
%   Inputs: M - Mach number
%   Outputs: D/q
%   taken from P.427 in Raymer textbook.
arguments (Input)
    M double
end

arguments (Output)
    D_q
end
Data = [0.525,0.3718
0.55,0.3737
0.575,0.3776
0.6,0.3839
0.625,0.3929
0.65,0.405
0.675,0.4206
0.7,0.4401
0.725,0.4634
0.75,0.4909
0.775,0.5225
0.8,0.5584
0.825,0.6017
0.85,0.6579
0.875,0.7303
0.9,0.8184
0.925,0.9243
0.95,1.0497
0.975,1.1973
1,1.3785
];

func =griddedInterpolant(Data(:,1),Data(:,2));

D_q = func(M);


end