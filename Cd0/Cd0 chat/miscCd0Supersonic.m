function Cd_misc = miscCd0Supersonic(M, params)
% Cd_misc = miscCd0Supersonic(M, params)
% Put your supersonic CD0 extras here (wave drag, store wave drag, etc.)

Cd_misc = 0;

% aim-120c    aim-9x      MK-83 JDAM

A_base = [19.63        19.63       153.9 ]; % using area of a circle0

d_q = 0.064 + (0.042*(M-3.84))*A_base;

% parasidic drag coefficient 

Cd_misc = Cd_misc + d_q / params.Sref; 

end
