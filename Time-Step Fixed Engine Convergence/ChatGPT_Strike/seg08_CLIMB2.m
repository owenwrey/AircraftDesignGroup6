function Tbl = seg08_CLIMB2(Tbl, istart, iend, p)
%SEG08_CLIMB2 Climb 1 to 30,000 ft at fixed Mach.

Tbl.Alt(istart:iend) = linspace(0,30000,(iend-istart+1));

for i = istart:iend
    [a, rho] = atmosAtAlt(Tbl.Alt(i), p);
    Tbl.rho(i) = rho;
    Tbl.TLapse(i) = p.TLapse(Tbl.Alt(i));

    Tbl.MACH(i) = 0.9;
    Tbl.KTAS(i) = Tbl.MACH(i) * a * p.mps2kts;
    Tbl.KEAS(i) = Tbl.KTAS(i) * sqrt(Tbl.rho(i)/p.rho_SL);

    Tbl.THROT(i)=1;
    Tbl.Thrust(i)=p.Thrust_SLS * Tbl.TLapse(i) * Tbl.THROT(i);
    Tbl.FF(i)=p.SFC_climb * Tbl.Thrust(i);

    % Aero must use i-1 weight fraction; ensure initialization
    Tbl = computeAeroFromState(Tbl, i, p);

    % Specific excess power / climb rate
    Tbl.dhdt(i) = ((Tbl.Thrust(i)*Tbl.KTAS(i)*p.kts2fps*60) - (Tbl.Drag(i)*Tbl.KTAS(i)*p.kts2fps*60)) / Tbl.Weight(i-1);
    Tbl.Ps(i) = Tbl.dhdt(i);

    Tbl.dTime(i) = (Tbl.Alt(i) - Tbl.Alt(i-1)) / Tbl.dhdt(i);
    Tbl.Time(i) = Tbl.Time(i-1) + Tbl.dTime(i);

    Tbl.FPA(i) = asind(Tbl.dhdt(i) / (Tbl.KTAS(i)*p.kts2fps*60));
    Tbl.GS(i)  = Tbl.KTAS(i) * cosd(Tbl.FPA(i));

    % distance credit
    Tbl.Dist(i) = Tbl.Dist(i-1) + (Tbl.GS(i) * Tbl.dTime(i) / 60);
    Tbl.dDist(i) = Tbl.Dist(i) - Tbl.Dist(i-1);

    if i==istart
        Tbl.dVdt(i)=0;
    else
        Tbl.dVdt(i) = (Tbl.KTAS(i)*p.kts2fps - Tbl.KTAS(i-1)*p.kts2fps) / (Tbl.dTime(i)*60);
    end

    Tbl.WtDrop(i)=0;
    Tbl = updateFuelAndWeight(Tbl, i, p);
end
end
