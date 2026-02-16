function Tbl = seg06_COMBAT(Tbl, istart, iend, p)
%SEG06_COMBAT Combat at SL, distance credit 100 NM, afterburner-ish profile.

Tbl.Dist(istart:iend) = linspace(Tbl.Dist(istart-1), Tbl.Dist(istart-1)+p.combat_nm, (iend-istart+1));

for i = istart:iend
    Tbl.dDist(i)=Tbl.Dist(i)-Tbl.Dist(i-1);

    Tbl.Alt(i)=0;
    [a, rho] = atmosAtAlt(0, p);
    Tbl.rho(i)=rho;
    Tbl.TLapse(i)=p.TLapse(0);

    Tbl.KTAS(i)=p.combat_KTAS;
    Tbl.KEAS(i)=Tbl.KTAS(i)*sqrt(Tbl.rho(i)/p.rho_SL);
    Tbl.MACH(i)=Tbl.KTAS(i)/(a*p.mps2kts);

    Tbl = computeAeroFromState(Tbl, i, p);

    % Original script solved for THROT via thrust=drag; keeping that behavior
    Tbl.Thrust(i)=Tbl.Drag(i);
    Tbl.THROT(i)=Tbl.Thrust(i)/(p.Thrust_SLS*Tbl.TLapse(i));
    Tbl.FF(i)=p.SFC_combat * Tbl.Thrust(i);

    Tbl.dhdt(i)=0; Tbl.Ps(i)=0;
    Tbl.FPA(i)=0;

    Tbl.dTime(i)=Tbl.dDist(i)/(Tbl.KTAS(i)/60);
    Tbl.Time(i)=Tbl.Time(i-1)+Tbl.dTime(i);

    Tbl.GS(i)=Tbl.KTAS(i);
    Tbl.dVdt(i)=0;

    Tbl.WtDrop(i)=0;
    Tbl = updateFuelAndWeight(Tbl, i, p);
end
end
