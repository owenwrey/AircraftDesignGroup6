function Cd_misc = miscCd0(M, params, regime)
% Cd_misc = miscCd0(M, params, regime)

switch regime
    case "subsonic"
        Cd_misc = miscCd0Subsonic(M, params);

    case "supersonic"
        Cd_misc = miscCd0Supersonic(M, params);

    otherwise
        error("Unknown regime: %s", regime);
end
end
