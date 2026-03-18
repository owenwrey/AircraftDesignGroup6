%% Weight Build Up the Sequel %%
clc; clear;

%% ================= GIVEN GEOMETRY =================
% ===== Fuselage =====
V_fuse = aircraft.fuselage.volume;     % Fuselage internal volume [ft^3]

% ===== Main Wing =====
S_w        = aircraft.wing.Area;  % Wing planform area [ft^2]
b_w        = aircraft.wing.span;   % Wingspan [ft]
c_root_w  = aircraft.wing.root;   % Root chord [ft]
c_tip_w   = aircraft.wing.tip;    % Tip chord [ft]
MAC_w     = aircraft.wing.MAC;   % Mean aerodynamic chord [ft]
Y_w       = aircraft.wing.MACloc;   % Spanwise MAC location [ft]

A         = b_w^2 / S_w;            % Wing aspect ratio [-]
lambda    = c_tip_w / c_root_w;     % Wing taper ratio [-]

% ===== Horizontal Tail =====
S_ht      = aircraft.ht.Area;  % Horizontal tail area [ft^2]
b_ht      = aircraft.ht.span;   % HT span [ft]
c_root_ht = aircraft.ht.;    % HT root chord [ft]
c_tip_ht  = aircraft.ht.;    % HT tip chord [ft]
MAC_ht    = aircraft.ht.;    % HT MAC [ft]
Y_ht      = aircraft.ht.;    % HT MAC location [ft]

% ===== Vertical Tail (Twin) =====
S_vt_total = aircraft.vt.Area; % Total vertical tail area (both tails) [ft^2]
S_vt       = aircraft.vt.span;  % Per-tail area [ft^2]
b_vt       = aircraft.vt.;   % Per-tail span [ft]
c_root_vt  = aircraft.vt.;  % VT root chord [ft]
c_tip_vt   = aircraft.vt.;   % VT tip chord [ft]
MAC_vt     = aircraft.vt.;   % VT MAC [ft]
Y_vt       = aircraft.vt.;   % VT MAC location [ft]

lambda_vt  = c_tip_vt / c_root_vt;  % VT taper ratio [-]

%% ================= GENERAL AIRCRAFT INPUTS =================
W_dg = (aircraft.);    % Design gross weight [lb]
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