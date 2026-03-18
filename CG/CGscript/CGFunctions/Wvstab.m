function [wvstab] = Wvstab(bv,Wto,nmax,Sv,Sr,MH,lv,Av,lamV,LamQuV)
%WVSTAB Estimates the weight of the vertical stablizer using method
%outlined by Roskam.
%   [wvstab] = Wvstab(bv,Wto,nmax,Sv,Sr,MH,lv,Av,lamV,LamQuV)
%   bv - vstab span
%   Wto - takeoff weight
%   nmax - max load factor
%   Sv - area vert stab
%   Sr - rudder area
%   MH - max mach at sea level
%   lv - distance from wing cbar/4 to vstab cbar/4
%   Av - aspect ratio v stab
%   lamV - labmda_v - taper ratio of vert stab
%   LamQuV - Lambda_1/4_V - quarter chord sweep of vert stab
%   angles in degrees.
arguments (Input)
    bv double
    Wto double
    nmax double
    Sv double
    Sr double
    MH double
    lv double
    Av double
    lamV double
    LamQuV double
end

arguments (Output)
    wvstab
end

zh = 0;

wvstab = 0.19*((1+zh/bv)^0.5 * (Wto*nmax)^0.363 * Sv^1.089 * MH^0.601 * ...
    lv^-0.726 * (1+Sr/Sv)^0.217 * Av^0.337 * (1 + lamV)^0.363 * cosd(LamQuV)^-0.484)^1.014;

end