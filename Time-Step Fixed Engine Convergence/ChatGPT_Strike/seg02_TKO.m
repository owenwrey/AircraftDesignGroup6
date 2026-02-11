function Tbl = seg02_TKO(Tbl, istart, iend, p)
%SEG02_TKO Takeoff (1 min). Afterburner (THROT = 1.25). No distance credit.

Tbl.Time(istart:iend) = linspace(Tbl.Time(istart-1), Tbl.Time(istart-1)+1, (iend-istart+1));

for i = istart:iend
    Tbl.dTime(i) = Tbl.Time(i) - Tbl.Time(i-1);

    Tbl.Alt(i) = 0;
    [~, rho] = atmosAtAlt(0, p);
    Tbl.rho(i) = rho;

    Tbl.KEAS(i)=0; Tbl.KTAS(i)=0; Tbl.MACH(i)=0;
    Tbl.FPA(i)=0; Tbl.GS(i)=0;
    Tbl.Dist(i)=Tbl.Dist(i-1); Tbl.dDist(i)=0;
    Tbl.dhdt(i)=0; Tbl.dVdt(i)=0;

    Tbl.CL(i)=0; Tbl.CD0(i)=0.025; Tbl.K1(i)=0; Tbl.K2(i)=0.05; Tbl.CDR(i)=0; Tbl.CD(i)=0;
    Tbl.L_D(i)=0; Tbl.Drag(i)=0;

    Tbl.TLapse(i)=1;
    Tbl.THROT(i)=1.25;
    Tbl.Thrust(i)=p.Thrust_SLS * Tbl.TLapse(i) * Tbl.THROT(i);
    Tbl.FF(i)=p.SFC_takeoff * Tbl.Thrust(i);

    Tbl.WtDrop(i)=0;
    Tbl = updateFuelAndWeight(Tbl, i, p);
end
end
