function aircraft = getK2_GI(aircraft, M, CL)
% WORK IN PROGRESS
%   Inputs: 
%       M - Mach number
%       CL - CL for current flight condition
%       CLdesign - wing design CL, from aircraft.constants.CLdesign
%   Outputs: K2 - K2 gridded interpolant
%
%   taken from P.446-449 in Raymer textbook.
% 

arguments (Input)
    aircraft struct
    M double
    CL double
end

arguments (Output)
    aircraft struct
end

CLdesign = aircraft.constants.CLdesign;


S = 0;

K0 = [
        0.2,0.2723
        0.25,0.2708
        0.3,0.2692
        0.35,0.2673
        0.4,0.2651
        0.45,0.2625
        0.5,0.2594
        0.55,0.2558
        0.6,0.2515
        0.65,0.2467
        0.7,0.2413
        0.75,0.2355
        0.8,0.2292
        0.85,0.2223
        0.9,0.2145
        0.95,0.2056
        1,0.1963
        1.05,0.189
        1.1,0.1856
        1.15,0.1864
        1.2,0.1915
        1.25,0.2011
        1.3,0.2137
        1.35,0.2275
        1.4,0.2418
        1.45,0.256
        1.5,0.2698
        1.55,0.2832
        1.6,0.297
        ];

K100 = [
        0.2,0.0922
        0.25,0.0923
        0.3,0.0921
        0.35,0.0919
        0.4,0.0918
        0.45,0.0917
        0.5,0.0919
        0.55,0.0923
        0.6,0.0931
        0.65,0.0936
        0.7,0.0933
        0.75,0.0919
        0.8,0.0905
        0.85,0.091
        0.9,0.0952
        0.95,0.1038
        1,0.1157
        1.05,0.1298
        1.1,0.1457
        1.15,0.1626
        1.2,0.1801
        1.25,0.1976
        1.3,0.214
        1.35,0.2283
        1.4,0.2418
        1.45,0.2554
        1.5,0.2692
        1.55,0.2831
        1.6,0.297
        ];

func0 = griddedInterpolant(K0(:,1),K0(:,2));
func100 = griddedInterpolant(K100(:,1),K100(:,2));

K2 = (S*(func100(M) + (1-S)*func0(M)));

aircraft.aero.K2_GI = K2;
end