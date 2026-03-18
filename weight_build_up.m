clc; clear;

%% ================= GIVEN GEOMETRY =================
% ===== Fuselage =====
V_fuse = 1145.11;     % Fuselage internal volume [ft^3]

% ===== Main Wing =====
S_w        = 741.75;  % Wing planform area [ft^2]
b_w        = 54.47;   % Wingspan [ft]
c_root_w  = 21.79;   % Root chord [ft]
c_tip_w   = 5.45;    % Tip chord [ft]
MAC_w     = 15.25;   % Mean aerodynamic chord [ft]
Y_w       = 10.89;   % Spanwise MAC location [ft]

A         = b_w^2 / S_w;            % Wing aspect ratio [-]
lambda    = c_tip_w / c_root_w;     % Wing taper ratio [-]

% ===== Horizontal Tail =====
S_ht      = 172.39;  % Horizontal tail area [ft^2]
b_ht      = 26.26;   % HT span [ft]
c_root_ht = 9.38;    % HT root chord [ft]
c_tip_ht  = 3.75;    % HT tip chord [ft]
MAC_ht    = 6.97;    % HT MAC [ft]
Y_ht      = 5.63;    % HT MAC location [ft]

% ===== Vertical Tail (Twin) =====
S_vt_total = 153.11; % Total vertical tail area (both tails) [ft^2]
S_vt       = 76.55;  % Per-tail area [ft^2]
b_vt       = 8.75;   % Per-tail span [ft]
c_root_vt  = 13.46;  % VT root chord [ft]
c_tip_vt   = 4.04;   % VT tip chord [ft]
MAC_vt     = 9.60;   % VT MAC [ft]
Y_vt       = 1.79;   % VT MAC location [ft]

lambda_vt  = c_tip_vt / c_root_vt;  % VT taper ratio [-]

%% ================= GENERAL AIRCRAFT INPUTS =================
W_dg = (0.6*34375);    % Design gross weight [lb]
N_z  = 8;     % Ultimate load factor [-]
M    = 1.2;      % Cruise Mach number [-]

%% ================= WING STRUCTURAL PARAMETERS =================
Lambda     = deg2rad(24);  % Wing quarter-chord sweep [rad]
t_c_root   = 0.12;         % Wing thickness-to-chord ratio at root [-]
S_csw      = 100;          % Wing control surface area [ft^2]

K_dw = 1.0;     % Wing dynamic pressure factor [-]
K_vs = 1.0;     % Wing variable sweep factor (fixed wing = 1) [-]

%% ================= EMPENNAGE GEOMETRY FACTORS =================
F_w  = 6.0;     % Fuselage width at HT intersection [ft]
B_h  = b_ht;    % HT span [ft]
H_t  = 5.0;     % HT height above fuselage centerline [ft]
H_v  = 10.0;    % VT height above fuselage centerline [ft]
L_t  = 15;      % Tail moment arm (wing AC to tail AC) [ft]

Lambda_vt = deg2rad(40);  % Vertical tail sweep angle [rad]
S_r  = 10;               % Rudder area [ft^2]
S_v  = S_vt;             % Vertical tail reference area [ft^2]
K_rht = 1.0;             % Rudder/HT correction factor [-]

%% ================= FUSELAGE GEOMETRY =================
K_dwf = 1.0;    % Fuselage design factor [-]
L_f   = 50;     % Fuselage length [ft]
D     = 6.0;    % Fuselage depth [ft]
W     = 6.0;    % Fuselage width [ft]

%% ================= LANDING GEAR =================
W_L = 0.95 * W_dg;  % Landing weight [lb]
N_l = 3;            % Ultimate landing load factor [-]
L_m = 4.5;          % Main gear length [ft]
L_n = 3.0;          % Nose gear length [ft]
N_w = 2;            % Number of wheels on nose gear [-]

K_cb  = 1.0;        % Gear structural factor [-]
K_tpg = 1.0;        % Gear type factor [-]

%% ================= ENGINE PARAMETERS =================
N_en = 2;           % Number of engines [-]
T    = 29000;       % Thrust per engine [lb]
T_e  = T;           % Equivalent thrust [lb]
W_en = 3980;        % Engine weight (each) [lb]

% ===== GE F110-GE-129 (spec-based) =====
D_e = 46.5/12;        % Engine maximum diameter [ft]  (46.5 in)  ≈ 3.875 ft

% Optional: overall engine length for sanity checks (not used in your equations)
L_engine = 181.9/12;  % Engine overall length [ft] (181.9 in) ≈ 15.16 ft

% ===== Installation-dependent (NOT engine-spec constants) =====
% These must come from your CAD / inlet + nacelle layout.
% If you don't have them, use placeholders but label them as assumptions.

L_d  = 8.0;   % Inlet/diffuser duct length: inlet lip -> engine face [ft] (airframe-dependent)
L_s  = 4.0;   % Accessory section length [ft] (modeling assumption unless you have packaging)
L_tp = 6.0;   % Tailpipe length: turbine exit -> nozzle exit [ft] (airframe-dependent)
L_sh = 15.2;  % Nacelle/shroud length covering the engine [ft] (should be ~ engine length)
L_ec = 12.0;  % Engine control routing length [ft] (airframe-dependent: cockpit/avionics -> engine)

K_vg = 1.0;   % Inlet/air induction installation factor [-] (assumption unless modeled)
K_d  = 1.0;   % Diffuser efficiency factor [-] (assumption unless modeled)
%% ================= SYSTEMS INPUTS =================
S_fw = 25;      % Firewall area [ft^2]
V_t  = 500;     % Total fuel volume [gal]
V_i  = 50;      % Integral tank volume [gal]
V_p  = 30;      % Self-sealing tank volume [gal]
N_t  = 4;       % Number of fuel tanks [-]
SFC  = 0.8;     % Specific fuel consumption [lb/(lbf·hr)]

S_cs = 200;     % Control surface area [ft^2]
N_s  = 4;       % Number of control systems [-]
N_c  = 1;     % Number of crew/passengers [-]
N_ci = 2;       % Number of crew interfaces [-]

K_vsh = 1.0;    % Hydraulic system factor [-]
N_u   = 3;      % Number of utility systems [-]

K_mc  = 1.0;    % Electrical system factor [-]
R_kva = 150;    % Electrical rating [kVA]
L_a   = 80;     % Electrical routing length [ft]
N_gen = 2;      % Number of generators [-]

W_uav = 800;    % Avionics weight estimate [lb]

%% ================= STRUCTURAL WEIGHTS =================
W_wing = 0.0103*K_dw*K_vs*(W_dg*N_z)^0.5 * ...
         S_w^0.622 * A^0.785 * t_c_root^(-1) * ...
         (1+lambda)^0.05 * (cos(Lambda))^(-1.0) * S_csw^0.04;

W_ht = 3.316*(1+(F_w/B_h))^(-2.0) * ((W_dg*N_z)/1000)^0.260 * S_ht^0.806;

W_vt = 0.452*K_rht*(1+H_t/H_v)^0.5 * (W_dg*N_z)^0.488 * ...
       S_vt^0.718 * M^0.341 * L_t^(-1.0) * ...
       (1+S_r/S_v)^0.348 * Lambda_vt^0.223 * ...
       (1+lambda_vt)^0.223 * (cos(Lambda_vt))^(-0.323);

W_fuselage = 0.499*K_dwf*W_dg^0.35*N_z^0.5 * ...
             L_f^0.25 * D^0.849 * W^0.685;

%% ================= LANDING GEAR WEIGHTS =================
W_main_gear = K_cb*K_tpg*(W_L/N_l)^0.25 * L_m^0.973;
W_nose_gear = (W_L/N_l)^0.290 * L_n^0.50 * N_w^0.525;

%% ================= ENGINE & SYSTEMS =================
W_engine_mounts = 0.013*N_en^0.795*T^0.579*N_z;

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
W_hydraulics  = 37.23*K_vsh*N_u^0.664;
W_electrical  = 172.2*K_mc*R_kva^0.152*N_c^0.10*L_a^0.10*N_gen^0.091;

W_avionics   = 2.117*W_uav^0.933;
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
W_fuel = 34375; 
W_gross = W_empty + W_fuel;
%% ================= DISPLAY =================

%% ================= AUTO-REPORT ALL WEIGHT VARIABLES =================
fprintf('\n===== ALL WEIGHT VARIABLES (W_*) (lb) =====\n')

% Get everything in the workspace
vars = whos;

% Loop through variables and print anything that starts with "W_"
for k = 1:length(vars)
    name = vars(k).name;

    if startsWith(name,'W_')
        val = eval(name);

        % Print scalars nicely
        if isnumeric(val) && isscalar(val)
            fprintf('%-20s = %12.2f\n', name, val);

        % If a vector shows up, print summary (shouldn't happen in this script)
        elseif isnumeric(val)
            fprintf('%-20s = [array size %s]\n', name, mat2str(size(val)));
        end
    end
end

fprintf('==========================================\n')

