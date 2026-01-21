function Cf = supersonicCf(Re, M)
% Cf = supersonicCf(Re, M)
% Paste your required supersonic turbulent skin-friction model here.

% PLACEHOLDER (same as subsonic right now):
Cf = 0.455 / ( (log10(Re))^2.58 * (1 + 0.144*M^2)^0.65 );

end
