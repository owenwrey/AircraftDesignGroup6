function aircraft = dimensionalize_aircraft(aircraft)
%DIMENSIONALIZE_AIRCRAFT  Compute aircraft dimensions from W0 and W/S.

%% Required checks 
req = @(s,f) isfield(s,f) && ~isempty(s.(f));

if ~isfield(aircraft,'constants') || ~req(aircraft.constants,'totalWeight') || ~req(aircraft.constants,'wingLoading')
    error('Need aircraft.constants.totalWeight and aircraft.constants.wingLoading.');
end
if ~isfield(aircraft,'fuselage') || ~req(aircraft.fuselage,'length') || ~req(aircraft.fuselage,'diameter')
    error('Need aircraft.fuselage.length and aircraft.fuselage.diameter.');
end
if ~isfield(aircraft,'wing') || ~req(aircraft.wing,'AR') || ~req(aircraft.wing,'taper_ratio')
    error('Need aircraft.wing.AR and aircraft.wing.taper_ratio.');
end
if ~isfield(aircraft,'ht') || ~req(aircraft.ht,'VolCoeff') || ~req(aircraft.ht,'AR') || ...
        ~req(aircraft.ht,'TaperRatio') || ~req(aircraft.ht,'leverArm_frac')
    error('Need aircraft.ht.VolCoeff, AR, TaperRatio, leverArm_frac.');
end
if ~isfield(aircraft,'vt') || ~req(aircraft.vt,'VolCoeff') || ~req(aircraft.vt,'AR') || ...
        ~req(aircraft.vt,'TaperRatio') || ~req(aircraft.vt,'leverArm_frac')
    error('Need aircraft.vt.VolCoeff, AR, TaperRatio, leverArm_frac.');
end

if ~isfield(aircraft.vt,'twinTail')
    aircraft.vt.twinTail = true;
end

%% Pull inputs
W0 = aircraft.constants.totalWeight;
WS = aircraft.constants.wingLoading;

Lf = aircraft.fuselage.length;
df = aircraft.fuselage.diameter;

ARw = aircraft.wing.AR;
TRw = aircraft.wing.taper_ratio;

CHT = aircraft.ht.VolCoeff;
ARh = aircraft.ht.AR;
TRh = aircraft.ht.TaperRatio;

CVT = aircraft.vt.VolCoeff;
ARv = aircraft.vt.AR;
TRv = aircraft.vt.TaperRatio;

Lh = aircraft.ht.leverArm_frac * Lf;
Lv = aircraft.vt.leverArm_frac * Lf;

%% Fuselage 
aircraft.fuselage.volume = pi*(df/2)^2 * Lf;

%% Wing 
S_w = W0 / WS;
b_w = sqrt(ARw * S_w);

Croot_w = (2 * S_w) / (b_w * (1 + TRw));
Ctip_w  = TRw * Croot_w;

MAC_w = ((2/3) * Croot_w * (1 + TRw + TRw^2)) / (1 + TRw);

aircraft.wing.Area = S_w;
aircraft.wing.span = b_w;
aircraft.wing.MAC  = MAC_w;
aircraft.wing.chord.root = Croot_w;
aircraft.wing.chord.tip  = Ctip_w;

%% Horizontal Tail 
S_HT = (CHT * MAC_w * S_w) / Lh;
b_HT = sqrt(ARh * S_HT);

Croot_HT = (2 * S_HT) / (b_HT * (1 + TRh));
Ctip_HT  = TRh * Croot_HT;

MAC_HT = ((2/3) * Croot_HT * (1 + TRh + TRh^2)) / (1 + TRh);

aircraft.ht.leverArm = Lh;
aircraft.ht.Area     = S_HT;
aircraft.ht.span     = b_HT;
aircraft.ht.MAC      = MAC_HT;
aircraft.ht.chord.root = Croot_HT;
aircraft.ht.chord.tip  = Ctip_HT;

%% Vertical Tail 
S_VT_total = (CVT * b_w * S_w) / Lv;   

if aircraft.vt.twinTail
    S_VT_each = S_VT_total / 2;
else
    S_VT_each = S_VT_total;
end

b_VT_each = sqrt(ARv * S_VT_each);

Croot_VT_each = (2 * S_VT_each) / (b_VT_each * (1 + TRv));
Ctip_VT_each  = TRv * Croot_VT_each;

MAC_VT_each = ((2/3) * Croot_VT_each * (1 + TRv + TRv^2)) / (1 + TRv);

aircraft.vt.leverArm  = Lv;
aircraft.vt.Area      = S_VT_total;   
aircraft.vt.Area_each = S_VT_each;
aircraft.vt.span_each = b_VT_each;
aircraft.vt.MAC_each  = MAC_VT_each;
aircraft.vt.chord.root_each = Croot_VT_each;
aircraft.vt.chord.tip_each  = Ctip_VT_each;

end

%{
%% Engine & Systems
K_vg = 1.62; %For variable geomtry (otherwise = 1)
K_d = 2.75; % Sqaure Inlet Duct
L_s_L_d = 1;
K_vsh = 1.425; %Variable sweep, 1 otherwise
K_mc = 1.45; %if mission completion required after failure, otherwise = 1

W_engine_mounts = 0.013*N_en^0.795*T^0.579*N_z;

W_firewall = 1.13*S_fw;
W_engine_section = 0.001*W_en^0.717*N_en*N_z;

W_air_induction = 13.29*K_vg*L_d^0.643*K_d^0.182 * N_en^1.498*(L_s_L_d)^(-0.373)*D_e;

W_tailpipe = 3.5*D_e*L_tp*N_en;
W_engine_cooling = 4.55*D_e*L_sh*N_en;
W_oil_cooling = 37.82*N_en^1.023;

W_engine_controls = 10.5*N_en^1.008*L_ec^0.222;
W_starter = 0.025*T_e^0.760*N_en^0.72;

W_fuel_system = 7.45*V_t^0.47 * (1+V_i/V_t)^(-0.095) * (1+V_p/V_t) * N_t^0.066 * N_en^0.052 * (T*SFC/1000)^0.249;

W_flight_controls = 36.28*M^0.003*S_cs^0.489*N_s^0.484*N_c^0.127;

W_instruments = 8.0 + 36.37*N_en^0.676*N_t^0.237 + 26.4*(1+N_ci)^1.356;
W_hydraulics  = 37.23*K_vsh*N_u^0.664;
W_electrical  = 172.2*K_mc*R_kva^0.152*N_c^0.10*L_a^0.10*N_gen^0.091;

W_avionics   = 2.117*W_uav^0.933;
W_furnishings = 217.6*N_c;
W_AC_AI = 3.2e-4*W_dg;

aircraft.enginesystems.numberOfEngines  = N_en;
aircraft.engine.thrust = T;   
aircraft.constants.n_ult = N_z; %1.5*limit load factor
aircraft.enginesystems.fwSurfaceArea = S_fw; %Firewall surface area
aircraft.enginesystems.engineWeight  = W_en;
aircraft.enginesystems.ductLength = L_d;
aircraft.enginesystems.engineDiameter = D_e;
aircraft.enginesystems.tailpipe = L_tp;
aircraft.enginesystems.shroudLength = L_sh;
aircraft.enginesystems.efcpDistance = L_ec; %engine front to cockpit distance (total if multiple engines)
aircraft.enginesystems.thrustPerEngine = T_e;
aircraft.enginesystems.totalFuelVolume = V_t; %gallons
aircraft.enginesystems.integralTankVolume = V_i; %gallons
aircraft.enginesystems.selfSealTankVolume = V_p; %gallons
aircraft.enginesystems.numberOfFuelTanks = N_t; 
aircraft.enginesystems.specficFuelConsumption = SFC; 
aircraft.enginesystems.mach = M;
aircraft.enginesystems.constrolSurfaceArea = S_cs; 
aircraft.enginesystems.numberOfFlightControlSystem = N_s; 
aircraft.enginesystems.numberOfCrew = N_c; 
aircraft.enginesystems.numberOfCrewEquivalence = N_ci; % 1 if one pilot, 1.2 if pilot + back seater, 2 if pilot and copilot
aircraft.enginesystems.numberOfHydraulicUtilityFunction = N_u; %Typically 5 to 15
aircraft.enginesystems.systemElectricalRating = R_kva; %Typically 110-160 for fighters and bombers
aircraft.enginesystems.electricalRoutingDistance =  L_a; %generators to avionics to cockpit, ft duct length, ft
aircraft.enginesystems.uninstalledAvionicsWeight = W_uav; %Typically 800-1400 lbs
aircraft.constants.totalWeight = W_dg; %typically 50-60% of internal fuel for military aircraft
%}