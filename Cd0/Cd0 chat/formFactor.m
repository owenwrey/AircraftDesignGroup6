function FF = formFactor(c, M, regime)
% FF = formFactor(c, M, regime)

% If FF is pre-defined (scalar), use it (fuselage/nacelle)
if isfield(c,'ff_scalar') && ~isempty(c.ff_scalar)
    FF = c.ff_scalar;
    return
end

% Wing/tail style FF from your formula
sweep = deg2rad(c.sweep_deg);
base  = (1 + (0.6/c.x_c)*c.t_c + 100*(c.t_c)^4);

% If you need a different supersonic FF, change this branch.
switch regime
    case "subsonic"
        corr = 1.34 * M^0.18 * cos(sweep)^0.28;
    case "supersonic"
        corr = 1.34 * M^0.18 * cos(sweep)^0.28; % placeholder
end

FF = base * corr;
end
