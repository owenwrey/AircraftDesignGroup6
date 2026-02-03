clc; clear;

%% ================= USER INPUTS =================
% -------- General ----------
W_dg   = 45000;     % Design gross weight [lb]
N_z    = 3.75;      % Ultimate load factor
M      = 0.9;       % Mach number

% -------- Wing -------------
S_w        = 600;    % Wing area [ft^2]
A          = 8.5;    % Aspect ratio
lambda     = 0.35;   % Taper ratio
Lambda     = deg2rad(35);  % Wing sweep [rad]
t_c_root   = 0.12;   % Thickness-to-chord at root
S_csw      = 100;    % Control surface area [ft^2]

K_dw = 1.0;
K_vs = 1.0;

% -------- Empennage --------
S_ht = 150;
S_vt = 120;
F_w  = 6.0;
B_h  = 20;
H_t  = 5.0;
H_v  = 10.0;
L_t  = 15;
Lambda_vt = deg2rad(40);
S_r  = 10;
S_v  = 120;
K_rht = 1.0;

% -------- Fuselage ---------
K_dwf = 1.0;
L_f   = 50;
D     = 6.0;
W     = 6.0;

% -------- Landing Gear -----
W_L = 0.95 * W_dg;
N_l = 3;
L_m = 4.5;
L_n = 3.0;
N_w = 2;

K_cb  = 1.0;
K_tpg = 1.0;

% -------- Engine -----------
N_en = 2;
T    = 18000;     % Thrust per engine [lb]
T_e  = T;
W_en = 2800;      % Engine weight [lb]

% -------- Engine Geometry --
D_e  = 4.0;
L_d  = 8.0;
L_s  = 4.0;
L_tp = 6.0;
L_sh = 5.0;
L_ec = 12.0;

K_vg = 1.0;
K_d  = 1.0;

% -------- Systems ----------
S_fw  = 25;
V_t   = 500;
V_i   = 50;
V_p   = 30;
N_t   = 4;
SFC   = 0.8;

S_cs = 200;
N_s  = 4;
N_c  = 120;
N_ci = 2;

K_vsh = 1.0;
N_u   = 3;

K_mc  = 1.0;
R_kva = 150;
L_a   = 80;
N_gen = 2;

W_uav = 800;

%% ================= STRUCTURE WEIGHTS =================

W_wing = 0.0103*K_dw*K_vs*(W_dg*N_z)^0.5 * ...
         S_w^0.622 * A^0.785 * t_c_root^(-1) * ...
         (1+lambda)^0.05 * (cos(Lambda))^(-1.0) * S_csw^0.04;

W_ht = 3.316*(1+F_w/B_h)^(-2.0) * ((W_dg*N_z)/1000)^0.260 * S_ht^0.806;

W_vt = 0.452*K_rht*(1+H_t/H_v)^0.5 * (W_dg*N_z)^0.488 * ...
       S_vt^0.718 * M^0.341 * L_t^(-1.0) * ...
       (1+S_r/S_v)^0.348 * Lambda_vt^0.223 * ...
       (1+lambda)^0.223 * (cos(Lambda_vt))^(-0.323);

W_fuselage = 0.499*K_dwf*W_dg^0.35*N_z^0.5 * L_f^0.25 * D^0.849 * W^0.685;

%% ================= LANDING GEAR =================

W_main_gear = K_cb*K_tpg*(W_L/N_l)^0.25 * L_m^0.973;
W_nose_gear = (W_L/N_l)^0.290 * L_n^0.50 * N_w^0.525;

%% ================= ENGINE & MOUNTS =================

W_engine_mounts = 0.013*N_en^0.795*T^0.579*N_z;

%% ================= SYSTEMS =================

W_firewall = 1.13*S_fw;

W_engine_section = 0.001*W_en^0.717*N_en*N_z;

W_air_induction = 13.29*K_vg*L_d^0.643*K_d^0.182 * ...
                  N_en^1.498*(L_s/L_d)^(-0.373)*D_e;

W_tailpipe = 3.5*D_e*L_tp*N_en;
W_engine_cooling = 4.55*D_e*L_sh*N_en;
W_oil_cooling = 37.82*N_en^1.023;

W_engine_controls = 10.5*N_en^1.008*L_ec^0.222;
W_starter = 0.025*T_e^0.760*N_en^0.72;

W_fuel_system = 7.45*V_t^0.47 * (1+V_i/V_t)^(-0.095) * ...
                (1+V_p/V_t) * N_t^0.066 * N_en^0.052 * ...
                (T*SFC/1000)^0.249;

W_flight_controls = 36.28*M^0.003*S_cs^0.489*N_s^0.484*N_c^0.127;

W_instruments = 8.0 + 36.37*N_en^0.676*N_t^0.237 + 26.4*(1+N_ci)^1.356;

W_hydraulics = 37.23*K_vsh*N_u^0.664;

W_electrical = 172.2*K_mc*R_kva^0.152*N_c^0.10*L_a^0.10*N_gen^0.091;

W_avionics = 2.117*W_uav^0.933;

W_furnishings = 217.6*N_c;

W_AC_AI = 3.2e-4*W_dg;

%% ================= TOTALS =================

W_structural = W_wing + W_ht + W_vt + W_fuselage + ...
               W_main_gear + W_nose_gear + W_engine_mounts;

W_systems = W_firewall + W_engine_section + W_air_induction + ...
            W_tailpipe + W_engine_cooling + W_oil_cooling + ...
            W_engine_controls + W_starter + W_fuel_system + ...
            W_flight_controls + W_instruments + W_hydraulics + ...
            W_electrical + W_avionics + W_furnishings + W_AC_AI;

W_empty = W_structural + W_systems;

%% ================= DISPLAY =================

fprintf('\nEMPTY WEIGHT BREAKDOWN (lb)\n')
fprintf('Wing: %.1f\n',W_wing)
fprintf('Fuselage: %.1f\n',W_fuselage)
fprintf('Empennage: %.1f\n',W_ht + W_vt)
fprintf('Landing Gear: %.1f\n',W_main_gear + W_nose_gear)
fprintf('Systems Total: %.1f\n',W_systems)
fprintf('----------------------------------\n')
fprintf('EMPTY WEIGHT TOTAL: %.1f lb\n',W_empty)
