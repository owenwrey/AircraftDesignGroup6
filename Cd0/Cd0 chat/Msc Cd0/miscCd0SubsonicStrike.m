function Cd_misc_Strike = miscCd0SubsonicStrike(M, params)
% Cd_misc = miscCd0Subsonic(M, params)

% aim-120c    aim-9x      MK-83 JDAM

A_base = [19.63        19.63       153.9 ]; % using area of a circle0

%d_q(1) = AIM120C
%d_q(2) = AIM-9x
d_q = 0.129+(0.419*(M-0.161)^2).*A_base;

% parasidic drag coefficient 

Cd_misc_Strike  = (4*d_q(3)+ 2*d_q(2)) / params.Sref ;

% call 2000lb bomb cluster on wing
% aim 9 missile and pylon

% 4 MK-83 JDAM 
% 2 AIm-9x
end