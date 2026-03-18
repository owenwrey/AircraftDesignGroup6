function [alpha] = TLapse(Alt,M,TR,selectEng)
%TLAPSE Thrust Lapse Calculation Based On Mattingly sec. 2.3.2
%   Calculates thrust lapse (alpha) from Altitude(ft.), Mach, Throttle Ratio,
%   and engine type. Taken from Mattingly pg. 38.
%
%   TR is 1 for engines running at MIL power (max no afterburner)
%   
%   Engine types:
%   HBPR - high bypasss Turbofan
%   LBPR - low bypass Turbofan
%   TURBOJET - turbojet
%
%   TLAPSE(Alt,M,TR,selectEng) = alpha
arguments (Input)
    Alt double 
    M double
    TR double
    selectEng string
end

arguments (Output)
    alpha double
end

% Neccessary values
[Ts,~,Ps,rhos,~,~] = atmosisa(0);
[Ts,~,Ps,~] = deal(Ts.*1.8,as./0.3048,Ps.*0.0208854,rhos.*0.0019403203);

[T,a,P,rho,~,~] = atmosisa(Alt.*0.3048);
[T,~,P,~] = deal(T.*1.8,a./0.3048,P.*0.0208854,rho.*0.0019403203);

th = T./Ts; th0 = th.*(1 + 0.2.*M.*M);
d = P./Ps; d0 = d.*(1 + 0.2.*M.*M).^(1.4./0.4);



% switch changes the equations used based on user engine input
switch selectEng

    case "HBPR" % High Bypass Ratio Turbofan

        if M >=0.9; error("Mach number too high for HBPR, see Mattingly p.38"); end;
        
        if (th0 <= TR)
            alpha = d0.*(1 - 0.49.*sqrt(M));
        else
            alpha = d0.*(1 - 0.49.*sqrt(M) - 3.*(th0 - TR)./(1.5 + M));
        end
        
    case "LBPR" % Low Bypass Ratio Turbofan
        
        if TR <= 1
            if (th0 <= TR)
                alpha = 0.6.*d0;
            else
                alpha = 0.6.*d0.*(1 - 3.8.*(th0-TR)./th0);
            end
        else % afterburning, TR > 1
            if (th0 <= TR)
                alpha = d0;
            else
                alpha = d0.*(1 - 3.5.*(th0-TR)./th0);
            end
        end

    case "TURBOJET" % self explanatory
        
        if TR > 1 % afterburner
            if (th0 <= TR)
                alpha = d0.*(1 - 0.3.*(th0-1)-0.1.*sqrt(M));
            else
                alpha = d0.*(1 - 0.3.*(th0-1)-0.1.*sqrt(M) - 1.5.*(th0-TR)./th0);
            end
        else % no afterburner
            if (th0 <= TR)
                alpha = 0.8.*d0.*(1 - 0.16.*sqrt(M));
            else
                alpha = 0.8.*d0.*(1 - 0.16.*sqrt(M) - 24.*(th0-TR)./(th0.*(9+M)));
            end
        end

    otherwise

        error("Engine Select Variable, please choose from one of the available engines.")

end


end