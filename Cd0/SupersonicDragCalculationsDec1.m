%% cd, k1, k2, cl calculating 
% cd0 supersonic , component 

clear all; 
close all; 
clc ;

%% initial values
alt = [25000];           % ft

mu = 1.789e-5;           % dynamic viscosity [kg/(m·s)]
m_sup  = 1.0998;         % supersonic mach number (design point)
k_default = 0.50e-5;     % ft  % polished sheet metal
sref = 578.846;          % reference wing area (ft^2)
cd_misc = 0;             % miscellaneous drag (to be updated later)
cd_l_p = 0;              % leakage + protuberance drag

%% component based analysis 

% wing
comp(1).name = 'wing';
comp(1).l    = 6.128;        % ft MAC (from board)
comp(1).swet = 1180.26;      % ft^2 wetted area (from board)
comp(1).t_c  = 0.045;        % thickness ratio t/c
comp(1).x_c  = 0.30;         % x/c location of max thickness
comp(1).sweep_angle_deg = 24;    % degrees
comp(1).sweep_angle = deg2rad(comp(1).sweep_angle_deg);
comp(1).ff = 1;              % supersonic FF = 1
comp(1).q = 1.0;             % q = 1 for mid-wing
comp(1).k = k_default;

% fuselage 
comp(2).name = 'fuselage';
comp(2).l    = 44.742;       % ft (board)
comp(2).d    = 15.805;       % ft diameter (board)
comp(2).swet = 1474.52;      % ft^2 wetted area (board computed)
comp(2).f    = comp(2).l / comp(2).d; 
comp(2).ff   = 1;            % supersonic FF = 1
comp(2).q    = 1.0;          % fuselage q ~ 1
comp(2).k    = k_default;

% horizontal tail
comp(3).name = 'horizontal tail';
comp(3).l    = 5.785;        % ft MAC from board
comp(3).swet = 240.46;       % ft^2 wetted area from board
comp(3).t_c  = 0.04;
comp(3).x_c  = 0.30;
comp(3).sweep_angle_deg = 30;
comp(3).sweep_angle = deg2rad(comp(3).sweep_angle_deg);
comp(3).ff = 1;              % supersonic FF = 1
comp(3).q = 1.03;            % pg 429
comp(3).k = k_default;

% vertical tail LEFT 
comp(4).name = 'vertical tail L';
comp(4).l    = 7.966;        % ft MAC (board)
comp(4).swet = 116.10;       % ft^2 wetted per tail
comp(4).t_c  = 0.30;
comp(4).x_c  = 0.30;
comp(4).sweep_angle_deg = 35;
comp(4).sweep_angle = deg2rad(comp(4).sweep_angle_deg);
comp(4).ff = 1;              % supersonic FF = 1
comp(4).q = 1.08;            % pg 429
comp(4).k = k_default;

% vertical tail RIGHT 
comp(5).name = 'vertical tail R';
comp(5).l    = 7.966;        % ft MAC
comp(5).swet = 116.10;       % ft^2 wetted
comp(5).t_c  = 0.30;
comp(5).x_c  = 0.30;
comp(5).sweep_angle_deg = 35;
comp(5).sweep_angle = deg2rad(comp(5).sweep_angle_deg);
comp(5).ff = 1;              % supersonic FF = 1
comp(5).q = 1.08;
comp(5).k = k_default;

% strut
comp(6).name = 'strut';
comp(6).l    = 3;        % ft
comp(6).swet = 5;        % ft^2
comp(6).t_c  = 0.12;
comp(6).x_c  = 0.30;
comp(6).sweep_angle_deg = 5;
comp(6).sweep_angle = deg2rad(comp(6).sweep_angle_deg);
comp(6).ff = 1;          % supersonic FF = 1
comp(6).q = 1.3;         % pg 429
comp(6).k = k_default;

% pylon 
comp(7).name = 'pylon';
comp(7).l    = 4;        % ft
comp(7).swet = 8;        % ft^2
comp(7).t_c  = 0.12;
comp(7).x_c  = 0.30;
comp(7).sweep_angle_deg = 0;
comp(7).sweep_angle = deg2rad(comp(7).sweep_angle_deg);
comp(7).ff = 1;          % supersonic FF = 1
comp(7).q = 1.4;         % pg 429
comp(7).k = k_default;

% nacelle 
comp(8).name = 'nacelle';
comp(8).l    = 12;       % ft
comp(8).d    = 3.2;      % ft
comp(8).swet = 120;      % ft^2
comp(8).f    = comp(8).l / comp(8).d;
comp(8).ff   = 1;        % supersonic FF = 1
comp(8).q    = 1.5;      % pg 425 - nacelle on wing
comp(8).k    = k_default;


%% supersonic CD0

cd0_sup = zeros(length(alt), length(m_sup));

for j = 1:length(alt)

    [temp_j, a_j, p_j, rho_j] = atmosisa(alt(j)*0.3048);  % a_j [m/s], rho_j [kg/m^3]

    for k = 1:length(m_sup)

        M = m_sup(k);
        V = M * a_j;    % true airspeed [m/s]

        cd0_sum = 0;   % accumulator

        % ---- mach dependent misc drag (missiles, racks, etc.) ----
        CD_misc_M = 0;
        % CD_misc_M = CD_misc_M + D_qAim(M)/sref;   % Aim-9 pylon drag
        % CD_misc_M = CD_misc_M + D_qClusterRack(M)/sref; % rack drag
        % CD_misc_M = CD_misc_M + D_q2KBwing(M) / sref; 

        % ---- leakage + protuberance drag ----
        cd_l_p = 0; % will be set after summing skin friction

        % ---- component loop (supersonic skin friction only) ----
        for i = 1:length(comp)

            l   = comp(i).l;      % ft
            sw  = comp(i).swet;   % ft^2
            k_c = comp(i).k;

            % FF and Q removed for supersonic
            FF = 1;
            Q = 1;

            % convert length to meters for Reynolds number
            l_m = l * 0.3048;     % m

            % Reynolds number (consistent SI units)
            Re = rho_j * V * l_m / mu;

            % supersonic Cf (flat plate)
            Cf_sup = 0.455 ./ (log10(Re).^2.58 .* (1 + 0.144*M^2).^0.65);

            % add component skin-friction cd0 (normalize with Sref in ft^2)
            cd0_sum = cd0_sum + (Cf_sup * FF * Q * sw) / sref;
        end

        % leakage/protuberance drag ≈ 12% of skin friction
        cd_l_p = 0.12 * cd0_sum;

        % wave drag (Sears-Haack approx for fuselage)
        d_fus = 5.4;          % ft (max fuselage diameter)
        ell   = 44.742;          % ft (fuselage length)
        Amax  = pi*(d_fus/2)^2;  % ft^2
        Ewd   = 3;               % empirical factor for fighters

        CD_wave = (9*pi/2)*(Amax/ell^2)^2 * Ewd;

        % total supersonic cd0
        cd0_sup(j,k) = cd0_sum + CD_misc_M + cd_l_p + CD_wave;

    end
end

cd0_sup
