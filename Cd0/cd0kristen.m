%% cd0 

% cd0 subsonic , component 

clear all; 
close all; 
clc ;

%% initial values

mu = 1.789e-5;           % dynamic viscosity [kg/(m·s)] (standard air, used with atmosisa outputs)

%% Inputs 

alt = [ ...
    0       % Start/Warmup
    0       % Takeoff
    25000   % Climb
    25000   % Cruise Outbound
    25000   % Descend to sea level 
    0       % combat
    25000   % Climb Return
    25000   % Cruise Inbound
    25000   % Descend Return
    0       % Loiter
    0   ];  % Landing
% 
% % KTAS from your Excel sheet (0 where blank)
% V_KTAS = [ ...
%     0
%     0
%     300
%     662
%     893.16
%     562
%     1.1851e3
%     861.9084
%     565
%     1.0428e3
%     278
%     0 ];
% 
% alt_m = alt * 0.3048;         % ft → m
% V_fts = V_KTAS * 1.68781;     % knots → ft/s
% 
% a_fts = zeros(size(alt));
% 
% for i = 1:length(alt_m)
%     [~, a_m_s, ~, ~] = atmosisa(alt_m(i));  % speed of sound in m/s
%     a_fts(i) = a_m_s * 3.28084;             % convert to ft/s
% 
% end
% 
% m = V_fts ./ a_fts;            % Mach number vector (12 values, one per mission segment)
m = [0 0 0.8 1.1 0.8 0.85 0.8 0.94 0.8 0.42 0];
% init values

k_default = 0.50e-5;     % ft  % polished sheet metal
sref = 578.846;          % reference area (ft^2)
cd_misc = 0;             % miscellaneous drag (to be updated later)
cd_l_p = 0;              % landing gear + protuberance drag (if used)

%% component based buildup 

% wing
comp(1).name = 'wing';
comp(1).l    = 6.128;        % ft (updated MAC from board)
comp(1).swet = 1180.26;      % ft^2 wetted area from board geometry
comp(1).t_c  = 0.045;        % thickness ratio t/c 
comp(1).x_c  = 0.30;         % x/c location of max thickness
comp(1).sweep_angle_deg = 24;    % degrees  (F/A-18 style LE sweep)
comp(1).sweep_angle = deg2rad(comp(1).sweep_angle_deg);
comp(1).ff = (1 + (0.6/comp(1).x_c)*comp(1).t_c + 100*(comp(1).t_c)^4) .* ...
              (1.34*m.^0.18 .* cos(comp(1).sweep_angle).^0.28);
comp(1).q = 1.0;         % pg 425, the f18 has a midwing, so q = 1
comp(1).k = k_default;

% fuselage 
comp(2).name = 'fuselage';
comp(2).l    = 50;       % ft fuselage length (board)
comp(2).d    = 12;       % ft max diameter (board)
comp(2).swet = 1474.52;      % ft^2 computed wetted area from board dims
comp(2).f    = comp(2).l / comp(2).d; 
comp(2).ff   = 0.9 + 5/comp(2).f^1.5 + comp(2).f/400;
comp(2).q    = 1.0;      % pg 429, the fuselage usually has negligible q
comp(2).k    = k_default;

% horizontal tail
comp(3).name = 'horizontal tail';
comp(3).l    = 5.785;       % ft MAC from board
comp(3).swet = 240.46;      % ft^2 wetted area from board geometry
comp(3).t_c  = 0.04;
comp(3).x_c  = 0.30;
comp(3).sweep_angle_deg = 30;    % realistic HT sweep for fighter
comp(3).sweep_angle = deg2rad(comp(3).sweep_angle_deg);
comp(3).ff = (1 + (0.6/comp(3).x_c)*comp(3).t_c + 100*(comp(3).t_c)^4) .* ...
              (1.34*m.^0.18 .* cos(comp(3).sweep_angle).^0.28);
comp(3).q = 1.03;        % pg 429
comp(3).k = k_default;

% vertical tail LEFT   (twin-tail configuration)
comp(4).name = 'vertical tail L';
comp(4).l    = 7.966;      % ft MAC (from board VT geometry)
comp(4).swet = 116.10;     % ft^2 wetted area per tail (from board)
comp(4).t_c  = 0.30;
comp(4).x_c  = 0.30;
comp(4).sweep_angle_deg = 35;
comp(4).sweep_angle = deg2rad(comp(4).sweep_angle_deg);
comp(4).ff = (1 + (0.6/comp(4).x_c)*comp(4).t_c + 100*(comp(4).t_c)^4) .* ...
              (1.34*m.^0.18 .* cos(comp(4).sweep_angle).^0.28);
comp(4).q = 1.08;         % pg 429
comp(4).k = k_default;

% vertical tail RIGHT  (twin-tail configuration)
comp(5).name = 'vertical tail R';
comp(5).l    = 7.966;      % ft MAC
comp(5).swet = 116.10;     % ft^2 wetted area
comp(5).t_c  = 0.30;
comp(5).x_c  = 0.30;
comp(5).sweep_angle_deg = 35;
comp(5).sweep_angle = deg2rad(comp(5).sweep_angle_deg);
comp(5).ff = (1 + (0.6/comp(5).x_c)*comp(5).t_c + 100*(comp(5).t_c)^4) .* ...
              (1.34*m.^0.18 .* cos(comp(5).sweep_angle).^0.28);
comp(5).q = 1.08;         % pg 429
comp(5).k = k_default;

% strut
comp(6).name = 'strut';
comp(6).l    = 3;        % ft
comp(6).swet = 5;        % ft^2
comp(6).t_c  = 0.12;
comp(6).x_c  = 0.30;
comp(6).sweep_angle_deg = 5;
comp(6).sweep_angle = deg2rad(comp(6).sweep_angle_deg);
comp(6).ff = (1 + (0.6/comp(6).x_c)*comp(6).t_c + 100*(comp(6).t_c)^4) .* ...
              (1.34*m.^0.18 .* cos(comp(6).sweep_angle).^0.28);
comp(6).q = 1.3;         % pg 429, wing strut has < drag than pylon
comp(6).k = k_default;

% pylon 
comp(7).name = 'pylon';
comp(7).l    = 4;        % ft
comp(7).swet = 8;        % ft^2
comp(7).t_c  = 0.12;
comp(7).x_c  = 0.30;
comp(7).sweep_angle_deg = 0;
comp(7).sweep_angle = deg2rad(comp(7).sweep_angle_deg);
comp(7).ff = (1 + (0.6/comp(7).x_c)*comp(7).t_c + 100*(comp(7).t_c)^4) .* ...
              (1.34*m.^0.18 .* cos(comp(7).sweep_angle).^0.28);
comp(7).q = 1.4;         % pg 429
comp(7).k = k_default;

% nacelle 
comp(8).name = 'nacelle';
comp(8).l    = 12;       % ft
comp(8).d    = 3.2;      % ft
comp(8).swet = 120;      % ft^2
comp(8).f    = comp(8).l / comp(8).d;
comp(8).ff   = 1 + 0.35/comp(8).f;
comp(8).q    = 1.5;      % pg 425 - nacelle mounted directly on wing
comp(8).k    = k_default;

%% cd0 combined (subsonic + supersonic) per mission segment
cd0_total = zeros(length(alt),1);

for j = 1:length(alt)

    % Atmos at segment altitude
    [~, a_j, ~, rho_j] = atmosisa(alt(j)*0.3048);  % a_j [m/s], rho_j [kg/m^3]

    M = m(j);
    V = M * a_j;    % [m/s]

    cd0_sum   = 0;
    CD_misc_M = 0;
    CD_wave   = 0;

    % CD_misc
           % aim-120c    aim-9x      MK-83 JDAM

    A_base = [19.63        19.63       153.9 ]; 
   
    for i = 1:length(comp)

        l   = comp(i).l;       % ft
        sw  = comp(i).swet;    % ft^2
        l_m = l * 0.3048;      % m

        Re = rho_j * V * l_m / mu;

        if M > 1
            % supersonic 
            FF = 1;
            Q  = 1;

            Cf = 0.455 ./ (log10(Re).^2.58 .* (1 + 0.144*M^2).^0.65);

            cd0_i = (Cf * FF * Q * sw) / sref;
            d_q = 0.064 +(0.042*(M-3.84))*A_base;
            Cd_misc_air_to_air  = (6*d_q(1)+ 2*d_q(2)) / sref ;
            Cd_misc_strike  = (4*d_q(3)+ 2*d_q(2)) / sref ;
        else
           % subsonic
            q_c = comp(i).q;
            k_c = comp(i).k;

            % FF selection (scalar vs Mach-array)
            if length(comp(i).ff) == 1
                ff_i = comp(i).ff;
            else
                ff_i = comp(i).ff(j);
            end

            % transition / roughness cap like you had
            Re_cut = 38.21 * (l/k_c)^1.053;
            Re_eff = min(Re, Re_cut);

            Cf = 0.455 / ( log10(Re_eff)^2.58 * (1 + 0.144*M^2)^0.65 );

            cd0_i = (Cf * ff_i * q_c * sw) / sref;
            d_q = 0.129+(0.419*(M-0.161)^2).*A_base;
            Cd_misc_air_to_air  = (6*d_q(1)+ 2*d_q(2)) / sref ;
            Cd_misc_strike  = (4*d_q(3)+ 2*d_q(2)) / sref ;
        end

        cd0_sum = cd0_sum + cd0_i;

    end

    % leakage + protuberance drag = 12% of skin friction sum
    cd_l_p = 0.12 * cd0_sum;

    % wave drag only when supersonic
    if M > 1
        d_fus = 5.4;                 % ft
        ell   = 44.742;              % ft
        Amax  = pi*(d_fus/2)^2;      % ft^2
        Ewd   = 3;

        CD_wave = (9*pi/2) * (Amax/ell^2)^2 * Ewd;
    end

    cd0_total(j) = cd0_sum + CD_misc_M + cd_l_p + CD_wave

end
