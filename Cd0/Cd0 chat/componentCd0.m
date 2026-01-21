function Cd0_i = componentCd0(c, M, rho, V, params, regime)
% Cd0_i = componentCd0(c, M, rho, V, params, regime)
%
% c.l_ft      characteristic length [ft]
% c.swet_ft2  wetted area [ft^2]
% c.q         interference factor [-]
% c.k_ft      roughness [ft]
% c.ff        form factor (scalar OR computed inside formFactor()) [-]

% Reynolds number
L_m = c.l_ft * 0.3048;            % [m]
Re  = rho * V * L_m / params.mu;  % [-]

% Roughness cutoff
Re = applyRoughnessCutoff(Re, c.l_ft, c.k_ft);

% Skin friction (sub vs sup)
Cf = skinFrictionCf(Re, M, regime);

% Form factor (kept same structure; change if your supersonic FF differs)
FF = formFactor(c, M, regime);

% Component CD0 contribution
Cd0_i = (Cf * FF * c.q * c.swet_ft2) / params.Sref;
end
