function [alpha] = TLapse(Alt,M,TR)
%TLAPSE Thrust Lapse Calculation Based On Mattingly sec. 2.3.2
%   Calculates thrust lapse (alpha) from Altitude(ft.), Mach, Throttle Ratio,
%   and engine type. Taken from Mattingly pg. 38.
%
%   TR is 1 for engines running at MIL power (max no afterburner)
%
%   TLAPSE(Alt,M,TR) = alpha

arguments (Input)
    Alt double 
    M double
    TR double
end

arguments (Output)
    alpha double
end

% Neccessary values
[Ts,as,Ps,rhos,~,~] = atmosisa(0);
[Ts,as,Ps,rhos] = deal(Ts.*1.8,as./0.3048,Ps.*0.0208854,rhos.*0.0019403203);

[T,a,P,rho,~,~] = atmosisa(Alt.*0.3048);
[T,a,P,rho] = deal(T.*1.8,a./0.3048,P.*0.0208854,rho.*0.0019403203);

th = T/Ts; th0 = th*(1 + 0.2*M*M);
d = P/Ps; d0 = d*(1 + 0.2*M*M)^(1.4/0.4);

        
if TR <= 1
    if (th0 <= TR)
        alpha = 0.6*d0;
    else
        alpha = 0.6.*d0.*(1 - 3.8*(th0-TR)/th0);
    end
else % afterburning, TR > 1
    if (th0 <= TR)
        alpha = d0;
    else
        alpha = d0.*(1 - 3.5*(th0-TR)/th0);
    end
end


end