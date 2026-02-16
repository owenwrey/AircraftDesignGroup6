function Tbl = seg01_SWT(Tbl, istart, iend, p)
%SEG01_SWT Start, warmup, ground operations (3 min). No distance credit.

Tbl.Time(istart:iend) = linspace(0,3,(iend-istart+1));

for i = istart:iend
    if i == istart
        Tbl.dTime(i) = 0;
        Tbl.Weight(i) = p.W0;
        Tbl.WtFrac(i) = 1;
        Tbl.FuelRem(i) = p.FuelReq;
    else
        Tbl.dTime(i) = Tbl.Time(i) - Tbl.Time(i-1);
    end

    Tbl.Alt(i)  = 0;
    [~, rho] = atmosAtAlt(Tbl.Alt(i), p);
    Tbl.rho(i) = rho;

    Tbl.KEAS(i) = 0; Tbl.KTAS(i)=0; Tbl.MACH(i)=0;
    Tbl.FPA(i)=0; Tbl.GS(i)=0;
    Tbl.Dist(i) = (i==istart)*0 + (i>istart)*Tbl.Dist(i-1);
    Tbl.dDist(i) = 0;
    Tbl.dhdt(i)=0; Tbl.dVdt(i)=0;

    Tbl.CL(i)=0; Tbl.CD0(i)=0.025; Tbl.K1(i)=0; Tbl.K2(i)=0.05; Tbl.CDR(i)=0; Tbl.CD(i)=0;
    Tbl.L_D(i)=0; Tbl.Drag(i)=0;

    Tbl.TLapse(i) = 1;
    Tbl.THROT(i)  = 0;
    Tbl.Thrust(i) = p.Thrust_SLS;
    Tbl.FF(i)     = p.SFC_idle * Tbl.Thrust(i);

    if i>istart
        Tbl.WtDrop(i)=0;
        Tbl = updateFuelAndWeight(Tbl, i, p);
    else
        Tbl.dFuel(i)=0;
        Tbl.FuelBurn(i)=0;
    end
end
end
