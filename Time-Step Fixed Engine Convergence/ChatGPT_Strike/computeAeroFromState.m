function Tbl = computeAeroFromState(Tbl, i, p)
%COMPUTEAEROFROMSTATE Computes CL/CD/Drag given Alt, KTAS, Weight history.
% Uses the same polar as the original script.

Tbl.CL(i)  = (Tbl.WtFrac(i-1) * p.W_S) / (0.5 * Tbl.rho(i) * (Tbl.KTAS(i) * p.kts2fps)^2);
Tbl.CD0(i) = 0.025;
Tbl.K1(i)  = 0.0;
Tbl.K2(i)  = 0.05;
Tbl.CDR(i) = 0.0;

Tbl.CD(i)  = Tbl.CD0(i) + Tbl.K1(i)*Tbl.CL(i) + Tbl.K2(i)*Tbl.CL(i)^2 + Tbl.CDR(i);
Tbl.L_D(i) = Tbl.CL(i) / Tbl.CD(i);

% Wing area is inferred from W0/W_S in your original script.
Tbl.Drag(i) = Tbl.CD(i) * 0.5 * Tbl.rho(i) * (Tbl.KTAS(i) * p.kts2fps)^2 * (p.W0 / p.W_S);
end
