function Tbl = seg11_LTR2(Tbl, istart, iend, p)
%SEG11_LTR2 Loiter at sea level for 20 min. No distance credit.

Tbl.Time(istart:iend) = linspace(Tbl.Time(istart-1), Tbl.Time(istart-1)+p.loiter_min, (iend-istart+1));

for i = istart:iend
    if i==istart
        Tbl.dTime(i)=0;
    else
        Tbl.dTime(i)=Tbl.Time(i)-Tbl.Time(i-1);
    end

    Tbl.Dist(i)=Tbl.Dist(i-1);
    Tbl.dDist(i)=0;

    Tbl.Alt(i)=0;
    [a, rho] = atmosAtAlt(0, p);
    Tbl.rho(i)=rho;
    Tbl.TLapse(i)=p.TLapse(0);

    Tbl.KTAS(i)=p.loiter_KTAS;
    Tbl.KEAS(i)=Tbl.KTAS(i)*sqrt(Tbl.rho(i)/p.rho_SL);
    Tbl.MACH(i)=Tbl.KTAS(i)/(a*p.mps2kts);

    Tbl = computeAeroFromState(Tbl, i, p);

    Tbl.Thrust(i)=Tbl.Drag(i);
    Tbl.THROT(i)=Tbl.Thrust(i)/(p.Thrust_SLS*Tbl.TLapse(i));
    Tbl.FF(i)=p.SFC_loiter*Tbl.Thrust(i);

    Tbl.dhdt(i)=0; Tbl.Ps(i)=0;
    Tbl.FPA(i)=0; Tbl.GS(i)=Tbl.KTAS(i); Tbl.dVdt(i)=0;

    Tbl.WtDrop(i)=0;
    if i>istart
        Tbl = updateFuelAndWeight(Tbl, i, p);
    else
        Tbl.dFuel(i)=0;
        Tbl.FuelBurn(i)=Tbl.FuelBurn(i-1);
        Tbl.FuelRem(i)=Tbl.FuelRem(i-1);
        Tbl.Weight(i)=Tbl.Weight(i-1);
        Tbl.WtFrac(i)=1;
    end
end
end
