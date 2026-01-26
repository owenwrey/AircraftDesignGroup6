function Cd0 = cd0Subsonic(alt_ft, M, comp, params)
% Cd0 = cd0Subsonic(alt_ft, M, comp, params)

[~, a, ~, rho] = atmosAtFt(alt_ft);
V = M * a; % [m/s]

Cd0_sum = 0;
for i = 1:numel(comp)
    Cd0_sum = Cd0_sum + componentCd0(comp(i), M, rho, V, params, "subsonic");
end

Cd_misc = miscCd0(M, params, "subsonic");
Cd_lp   = 0.12 * Cd0_sum;

Cd0 = Cd0_sum + Cd_misc + Cd_lp;
end
