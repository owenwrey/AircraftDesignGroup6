function Tbl = seg05_DESC1(Tbl, istart, iend, p)
%SEG05_DESC1 Descent 1 at -3000 ft/min to sea level.

Tbl.Alt(istart:iend) = linspace(Tbl.Alt(istart-1), 0, (iend-istart+1));

for i = istart:iend
    [a, rho] = atmosAtAlt(Tbl.Alt(i), p);
    Tbl.rho(i)=rho;
    Tbl.TLapse(i)=p.TLapse(Tbl.Alt(i));

    Tbl.MACH(i)=0.9;
    Tbl.KTAS(i)=Tbl.MACH(i)*a*p.mps2kts;
    Tbl.KEAS(i)=Tbl.KTAS(i)*sqrt(Tbl.rho(i)/p.rho_SL);

    Tbl.THROT(i)=1;
    Tbl.Thrust(i)=p.Thrust_SLS*Tbl.TLapse(i)*Tbl.THROT(i);
    Tbl.FF(i)=p.SFC_descent*Tbl.Thrust(i);

    Tbl.dhdt(i)=p.descent_fpm;
    Tbl.Ps(i)=Tbl.dhdt(i);

    Tbl = computeAeroFromState(Tbl, i, p);

    Tbl.dTime(i) = (Tbl.Alt(i) - Tbl.Alt(i-1)) / Tbl.dhdt(i);
    Tbl.Time(i)  = Tbl.Time(i-1) + Tbl.dTime(i);

    Tbl.FPA(i) = asind(Tbl.dhdt(i) / (Tbl.KTAS(i)*p.kts2fps*60));
    Tbl.GS(i)  = Tbl.KTAS(i) * cosd(Tbl.FPA(i));

    Tbl.Dist(i)  = Tbl.Dist(i-1) + Tbl.GS(i)*Tbl.dTime(i)/60;
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
