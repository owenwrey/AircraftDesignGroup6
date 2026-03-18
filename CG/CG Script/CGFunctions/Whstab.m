function [Whorizontal] = Whstab(Wto,nmax,Sh,bh,trh,cbar,lh)
%Whstab horizontal stablizer weight estimation
%   [Whorizontal] = Whstab(Wto,nmax,Sh,bh,trh,cbar,lh)
%   calculates weight of horizontal stblizer using GD method (Roksam)
%   Wto - weight takeoff
%   nmax - max load factor
%   Sh - area of horiz. stab.
%   bh - span of horiz. stab.
%   trh - horizontal tail max root thickness
%   cbar - mean aerodynamic chord
%   lh - length from main wing cbar/4 to Hstab cbar/4

arguments (Input)
    Wto double
    nmax double
    Sh double
    bh double
    trh double
    cbar double
    lh double
end

arguments (Output)
    Whorizontal double
end

Whorizontal = 0.0034*( (Wto*nmax)^0.813 * Sh^0.584 * (bh/trh)^0.033 * (cbar/lh)^0.28)^0.915;

end