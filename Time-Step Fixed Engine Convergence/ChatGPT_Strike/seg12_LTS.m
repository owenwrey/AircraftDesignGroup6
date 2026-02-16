function Tbl = seg12_LTS(Tbl, istart, iend, p)
%SEG12_LTS Landing, ground operations, shutdown (3 min). No distance credit.

Tbl.Time(istart:iend) = linspace(Tbl.Time(istart-1), Tbl.Time(istart-1)+3, (iend-istart+1));

for i = istart:iend
    if i==istart
        Tbl.dTime(i)=0;
    else
        Tbl.dTime(i)=Tbl.Time(i)-Tbl.Time(i-1);
    end

    Tbl.Alt(i)=0;
    [~, rho] = atmosAtAlt(0, p);
    Tbl.rho(i)=rho;

    Tbl.KEAS(i)=0; Tbl.KTAS(i)=0; Tbl.MACH(i)=0;
    Tbl.FPA(i)=0; Tbl.GS(i)=0;
    Tbl.Dist(i)=Tbl.Dist(i-1); Tbl.dDist(i)=0;
    Tbl.dhdt(i)=0; Tbl.dVdt(i)=0;

    Tbl.CL(i)=0; Tbl.CD0(i)=0.025; Tbl.K1(i)=0; Tbl.K2(i)=0.05; Tbl.CDR(i)=0; Tbl.CD(i)=0;
    Tbl.L_D(i)=0; Tbl.Drag(i)=0;

    Tbl.TLapse(i)=1;
    Tbl.THROT(i)=1;
    Tbl.Thrust(i)=p.Thrust_SLS*Tbl.TLapse(i)*Tbl.THROT(i);
    Tbl.FF(i)=p.SFC_shutdown*Tbl.Thrust(i);

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
