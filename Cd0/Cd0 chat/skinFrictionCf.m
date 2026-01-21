function Cf = skinFrictionCf(Re, M, regime)
% Cf = skinFrictionCf(Re, M, regime)

switch regime
    case "subsonic"
        Cf = 0.455 / ( (log10(Re))^2.58 * (1 + 0.144*M^2)^0.65 );

    case "supersonic"
        Cf = supersonicCf(Re, M); % YOU replace inside supersonicCf.m

    otherwise
        error("Unknown regime: %s", regime);
end
end
