function Cd_misc = miscCd0Supersonic(M, params)
% Cd_misc = miscCd0Supersonic(M, params)
% Put your supersonic CD0 extras here (wave drag, store wave drag, etc.)

Cd_misc = 0;

% Keep your existing term if it still applies:
Cd_misc = Cd_misc + D_q2KBwing(M) / params.Sref;

% Add wave drag here if your method has it:
% Cd_misc = Cd_misc + Cd_wave(M, geom, ...);

end
