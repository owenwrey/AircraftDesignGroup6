function aircraft = EWB(aircraft)

% initialize variables
W0 = aircraft.weight.total;
T2C = aircraft.wing.T2C;                    % thickness to chord ratio
Nz = aircraft.constants.ultLoad;            % ultimate load factor
M = aircraft.constants.maxMach;             % max Mach number

% wing
wing.area = aircraft.wing.Area;
wing.AR = aircraft.wing.AR;
wing.QuarterChordSweep = aircraft.wing.sweep;
wing.TaperRatio = aircraft.wing.taper_ratio;
wing.span = aircraft.wing.span;
rootChord = aircraft.wing.chord.root;
tipChord = aircraft.wing.chord.tip;
chordFrac = .2;                             % control surface chord fraction

% calculate wing control surface areas
% ailerons: 20% local chord, going from 50% to 90% span
% flaps: 20% local chord, going from 0% to 50% span
chord50 = (rootChord-tipChord)*.5 + tipChord;       % chord @ 50% span
chord90 = (rootChord-tipChord)*.9 + tipChord;       % chord @ 90% span
chord0_20 = rootChord*chordFrac;                    % 20% chord @ root
chord50_20 = chord50*chordFrac;                     % 20% chord @ 50% span
chord90_20 = chord90*chordFrac;                     % 20% chord @ 90% span
span50 = wing.span*.5;                              % 50% wingspan
span40 = wing.span*.4;                              % 40% wingspan
aircraft.wing.flapArea = .5*(chord0_20 + chord50_20)*span50;        % flap area, ft^2
aircraft.wing.aileronArea = .5*(chord90_20 + chord50_20)*span40;    % aileron area, ft^2
S_csw = aircraft.wing.flapArea + aircraft.wing.aileronArea;         % total control surface area

aircraft.wing.controlSurfaceArea = S_csw;

% horizontal tail
ht.area = aircraft.ht.Area;
ht.QuarterChordSweep = aircraft.ht.sweep;
ht.AspectRatio = aircraft.ht.AR;
ht.TaperRatio = aircraft.ht.TaperRatio;

% vertical tail
vt.area = aircraft.vt.Area;
vt.QuarterChordSweep = aircraft.vt.sweep;
vt.AspectRatio = aircraft.vt.AR;
vt.TaperRatio = aircraft.vt.TaperRatio;

% 25% chord and 90% span
% calculate rudder area, 25% local chord, going from 0% to 90% span
vt.rootChord = 1;
vt.chord90 = 1;
vt.chord0_20 = vt.rootChord*.2;
vt.chord90_20 = vt.chord90*.2;
vt.span = sqrt(vt.AspectRatio*vt.area);
aircraft.vt.rudderArea = .5*(vt.chord0_20 + vt.chord90_20)*vt.span*.9;


%% Wing, ht, vt, fuselage and gear weights
% from Raymer Fighter/Attack Weights, Sec 15.3.1

% main wing
aircraft.wing.weight = .0103*(W0*Nz)^(.5) * (wing.area)^(.622) * (wing.AR)^(.785)...
                        * T2C^(-.4) * (1+wing.TaperRatio)^(.05) * ...
                        (cosd(wing.QuarterChordSweep))^(-1) * (S_csw)^(.04);

% horizontal tail
Fw = 3.2;   % fuselage width at ht intersection
Bh = sqrt(wing.AR*wing.area);

aircraft.ht.weight = 3.316*(1+Fw/Bh)^(-2) * (W0*Nz/1000)^(.26) * ht.area^(.806);

% vertical tail
K_rht = 1.047;
L_t = aircraft.vt.leverArm;

aircraft.vt.weight = .452*K_rht * (W0*Nz)^(.488) * vt.area^(.718) * M^(.341)...
    * L_t^(-1) * (1+aircraft.vt.rudderArea/vt.area)^(.348) * vt.AspectRatio^.223...
    * (1+vt.TaperRatio)^(.25) * cosd(vt.QuarterChordSweep)^(-.323);

% fuselage
K_dwf = 1;                                  % delta wing multiplier
L = aircraft.fuselage.length;
D = 6;      % fuselage structural depth, ft
W = 6;      % fuselage structural width, ft

aircraft.fuselage.weight = .499*K_dwf * W0^(.35) * Nz^(.25) * L^(.5) * D^(.849)...
    * W^(.685);

% main landing gear
N_gear = 5.5;
K_cb = 1;                                   % cross-beam gear (we are not cross-beam)
K_tpg = .826;                               % tripod gear (we are tripod)
W_l = aircraft.weight.totalOnLanding;       % landing design gross weight, lb
N_l = 1.5*N_gear;                           % ultimate landing load factor
L_m = aircraft.gear.mg.extendedLength;      % extended length of main gear, in

aircraft.gear.mg.weight = K_cb * K_tpg * (W_l*N_l)^(.25) * L_m^(.973);

% nose landing gear
N_nw = 2;                                   % number of nosewheels
L_n = aircraft.gear.ng.extendedLength;      % extended length of nose gear, in

aircraft.gear.ng.weight = (W_l*N_l)^(.29) * L_n^(.5) * N_nw^(.525);




%% Engine & Systems
K_vg = 1.62;            % For variable geomtry (otherwise = 1)
K_d = 2.75;             % Sqaure Inlet Duct
L_s_L_d = 1;            % single duct length to duct length, 1 bc each engine has its own duct
L_d = 8;                % duct length
K_vsh = 1.425;          % Variable sweep, 1 otherwise
K_mc = 1.45;            % if mission completion required after failure, otherwise = 1
W_uav = 800;            % uninstalled avionics weight (800 lbs maybe)
S_fw = ((46.5/12)*(182/12))*3 + (47/12)^2-(pi*(46.5/12)^2)*2;       % firewall surface area, ft^2
S_cs = S_csw + ht.area + aircraft.vt.rudderArea;        % total control surface area

N_en = 2;               % number of engines
N_c = 1;                % number of crew
N_s = 2;                % number of flight control systems
N_ci = 1;               % number of crew equivalents
N_u = 9;                % # hydraulic functions (gear, ailerons, flaps, ht, vt, airbrake, brakes, canopy, folding wings)
N_gen = N_en;           % number of generators

W_en = aircraft.engine.weight;
T_e = aircraft.engine.thrust;               % thrust per engine
T = T_e * N_en;                             % total thrust
D_e = 46.5/12;          % engine diameter, ft
L_tp = 0;               % length of tailpipe, ft. 0 bc engine has nozzle???
L_sh = 3;               % length of cooling shroud, ft, just a guess
L_ec = (aircraft.engine.cg.x - aircraft.cockpit.cg.x) * N_en;     % dist from engine to cockpit, total if mult engines, ft
SFC = aircraft.engine.TSFC(1,1);
R_kva = 140;                                % system electrical rating, 110-160 for fighters
L_a = aircraft.engine.cg.x - aircraft.cockpit.cg.x;               % electrical routing dist, ft (between gen, avionics, & cockpit)

V_t = aircraft.constants.fuelVolume;    % fuel volume, gal
V_i = aircraft.fuelSys.VI;              % integral fuel volume, gal
% tanks made from cavities within the airframe, not dedicated fuel tanks
V_p = aircraft.fuelSys.VP;              % self-sealing tanks volume, gal
N_t = aircraft.fuelSys.Nt;

% engine structure calcs
engine.mounts = 0.013*N_en^(0.795) * T^(0.579) * Nz;

engine.firewall = 1.13*S_fw;

engine.section = 0.01*W_en^(0.717) * N_en * Nz;

engine.airInduction = 13.29 * K_vg * L_d^(0.643) * K_d^(0.182) * N_en^(1.498) ...
    * (L_s_L_d)^(-0.373) * D_e;

engine.tailpipe = 3.5 * D_e * L_tp * N_en;

engine.cooling = 4.55 * D_e * L_sh * N_en;

engine.oilCooling = 37.82 * N_en^(1.023);

engine.controls = 10.5 * N_en^(1.008) * L_ec^(0.222);

engine.starter = 0.025 * T_e^(0.760) * N_en^(0.72);

% misc system calcs
misc.fuelSystem = 7.45 * V_t^(0.47) * (1+V_i/V_t)^(-0.095) * (1+V_p/V_t)...
    * N_t^(0.066) * N_en^(0.052) * (T*SFC/1000)^(0.249);

misc.flightControls = 36.28 * M^(0.003) * S_cs^(0.489) * N_s^(0.484) * N_c^(0.127);

misc.instruments = 8.0 + 36.37 * N_en^(0.676) * N_t^(0.237) + 26.4*(1+N_ci)^(1.356);

misc.hydraulics  = 37.23 * K_vsh * N_u^(0.664);

misc.electrical  = 172.2 * K_mc * R_kva^(0.152) * N_c^(0.10) * L_a^(0.10) * ...
    N_gen^(0.091);

misc.avionics   = 2.117 * W_uav^(0.933);

misc.furnishings = 217.6 * N_c;

misc.AC_AI = 201.6 * ((W_uav+200*N_c)/1000)^(.735);

misc.handlingGear = 3.2e-4 * W0;


aircraft.engine.firewallSA = S_fw;
aircraft.aero.totalControlSurfaceArea = S_cs;
aircraft.enginesystems.numberOfEngines  = N_en;
aircraft.enginesystems.engineDiameter = D_e;
aircraft.enginesystems.tailpipe = L_tp;
aircraft.enginesystems.shroudLength = L_sh;
aircraft.enginesystems.efcpDistance = L_ec;         % engine front to cockpit distance (total if multiple engines)
aircraft.enginesystems.numberOfFlightControlSystem = N_s; 


% prevWeight = aircraft.weight.empty;

aircraft.weight.empty = aircraft.wing.weight + aircraft.ht.weight + aircraft.engine.weight*2 ... 
    + aircraft.engine.nozzle.weight*2 + aircraft.vt.weight + aircraft.fuselage.weight + aircraft.gear.mg.weight ...
    + aircraft.gear.ng.weight + sum(cell2mat(struct2cell(engine))) + sum(cell2mat(struct2cell(misc)));


aircraft.engine.structures = engine;
aircraft.weight.misc = misc;

% weightDiff = prevWeight - aircraft.weight.empty;

% weights = [weightDiff;aircraft.weight.empty];


end