function Tbl = seg09_CRUISE_IBD(Tbl, istart, iend, p)
%SEG09_CRUISE_IBD Cruise outbound: 30k ft, 650 NM.

Tbl.Dist(istart:iend) = linspace(Tbl.Dist(istart-1), Tbl.Dist(istart-1)+p.cruiseIn_nm, (iend-istart+1));

for i = istart:iend
    Tbl.dDist(i) = Tbl.Dist(i) - Tbl.Dist(i-1);

    Tbl.Alt(i) = 30000;
    [a, rho] = atmosAtAlt(Tbl.Alt(i), p);
    Tbl.rho(i)=rho;
    Tbl.TLapse(i)=p.TLapse(Tbl.Alt(i));

    Tbl.KTAS(i)=p.cruiseIn_KTAS;
    Tbl.KEAS(i)=Tbl.KTAS(i)*sqrt(Tbl.rho(i)/p.rho_SL);
    Tbl.MACH(i)=Tbl.KTAS(i)/(a*p.mps2kts);

    Tbl = computeAeroFromState(Tbl, i, p);

    % Level flight: thrust = drag
    Tbl.Thrust(i)=Tbl.Drag(i);
    Tbl.THROT(i)=Tbl.Thrust(i)/(p.Thrust_SLS*Tbl.TLapse(i));
    Tbl.FF(i)=p.SFC_cruise * Tbl.Thrust(i);

    Tbl.dhdt(i)=0; Tbl.Ps(i)=0;
    Tbl.FPA(i)=0;

    Tbl.dTime(i) = Tbl.dDist(i) / (Tbl.KTAS(i)/60);
    Tbl.Time(i)  = Tbl.Time(i-1) + Tbl.dTime(i);

    Tbl.GS(i)=Tbl.KTAS(i);
    Tbl.dVdt(i)=0;

    Tbl.WtDrop(i)=0;
    Tbl = updateFuelAndWeight(Tbl, i, p);
end
end
