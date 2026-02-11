function Tbl = updateFuelAndWeight(Tbl, i, p)
%UPDATEFUELANDWEIGHT Updates dFuel, FuelBurn, FuelRem, Weight, WtFrac.

Tbl.dFuel(i) = (Tbl.FF(i) * Tbl.dTime(i)) / 60;  % lb
Tbl.FuelBurn(i) = sum(Tbl.dFuel(1:i));
Tbl.FuelRem(i)  = Tbl.FuelRem(i-1) - Tbl.dFuel(i);
Tbl.Weight(i)   = Tbl.Weight(i-1) - Tbl.dFuel(i) - Tbl.WtDrop(i);

if Tbl.Weight(i-1) ~= 0
    Tbl.WtFrac(i) = Tbl.Weight(i) / Tbl.Weight(i-1);
else
    Tbl.WtFrac(i) = 1;
end
end
