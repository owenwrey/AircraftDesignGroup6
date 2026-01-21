function Cd_misc = miscCd0Subsonic(M, params)
% Cd_misc = miscCd0Subsonic(M, params)

Cd_misc = 0;

% Your current term:
Cd_misc = Cd_misc + D_q2KBwing(M) / params.Sref;

end
