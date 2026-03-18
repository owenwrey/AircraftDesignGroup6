function [Wfuselage] = Wfuse(Wto,lf,hf)
%WFUSE calculates weight of fuselage from roksam 5.28
%   qbarL - 
arguments (Input)
    Wto
    lf
    hf
end

arguments (Output)
    Wfuselage
end

kinl = 1.25; %because inlet is on fuselage
qbarL = 1.0536; % slugs/ft*s^2

Wfuselage = 11.03*(kinl)^1.23 * (qbarL/100)^0.245 * (Wto/1000)^0.98 * (lf/hf)^0.61;


end