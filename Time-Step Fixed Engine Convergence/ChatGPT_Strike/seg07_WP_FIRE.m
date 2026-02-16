function Tbl = seg07_WP_FIRE(Tbl, istart, iend, p)
%SEG07_WP_FIRE Weapons fire / drop (instantaneous point, no time or distance).
% Payload drop is modeled here.

for i = istart:iend
    Tbl.Time(i)  = Tbl.Time(i-1);
    Tbl.dTime(i) = 0;
    Tbl.Dist(i)  = Tbl.Dist(i-1);
    Tbl.dDist(i) = 0;

    Tbl.Alt(i)=0;
    [a, rho] = atmosAtAlt(0, p);
    Tbl.rho(i)=rho;
    Tbl.TLapse(i)=p.TLapse(0);

    Tbl.KTAS(i)=p.combat_KTAS;
    Tbl.KEAS(i)=Tbl.KTAS(i)*sqrt(Tbl.rho(i)/p.rho_SL);
    Tbl.MACH(i)=Tbl.KTAS(i)/(a*p.mps2kts);

    Tbl = computeAeroFromState(Tbl, i, p);

    Tbl.THROT(i)=1.25;
    Tbl.Thrust(i)=p.Thrust_SLS*Tbl.TLapse(i)*Tbl.THROT(i);
    Tbl.FF(i)=p.SFC_wtDrop * Tbl.Thrust(i);

    Tbl.dhdt(i)=0; Tbl.Ps(i)=0;
    Tbl.FPA(i)=0; Tbl.GS(i)=Tbl.KTAS(i); Tbl.dVdt(i)=0;

    % Drop stores here
    Tbl.WtDrop(i)=p.payloadDrop_lb;

    % With dTime=0, fuel burn is zero for this row (keep it that way)
    Tbl.dFuel(i)=0;
    Tbl.FuelBurn(i)=Tbl.FuelBurn(i-1);
    Tbl.FuelRem(i)=Tbl.FuelRem(i-1);

    Tbl.Weight(i)=Tbl.Weight(i-1) - Tbl.WtDrop(i);
    Tbl.WtFrac(i)=1;
end
end
