clear all; clc;

% aircraft = struct;
aircraft = struct;

% Geometry

% Empty Weight Buildup

% CG and Inertia
aircraft.cg.x = 20;
aircraft.cg.y = 0;
aircraft.cg.z = 10;

% Landing Gear

% Time-Step Mission
aircraft.weight.total = 60e3;
aircraft.weight.fuel = 20e3;
aircraft.weight.empty = 40e3;
aircraft.weight.totalOnLanding = 42e3;

%% -| Aircraft Struct (constant "variables")|------------------------------
% All variables that do not change in the loop but are necessary for
% calculations in the loop (e.g. aircraft.constants.wingLoading, aircraft.weight.tolerance)

% General
aircraft.constants.wingLoading = 80; % [lbf/ft]
aircraft.engine.weight = 5000;
aircraft.engine.thrust = 35000;
aircraft.engine.thrustMil = 26000;
aircraft.engine.TSFC = .67;

aircraft.constants.fuelVolume = 3500;
aircraft.fuelSys.VI = 0;                % integral fuel volume, gal
aircraft.fuelSys.VP = aircraft.constants.fuelVolume/2;  % self-sealing tanks volume, gal
aircraft.fuelSys.Nt = 4;

% Geometry
aircraft.fuselage.length   = 48;
aircraft.fuselage.diameter = 6;
aircraft.fuselage.cg.x = 24;
aircraft.fuselage.cg.y = 0;
aircraft.fuselage.cg.z = 10;

aircraft.wing.AR          = 3;
aircraft.wing.taper_ratio = 0.2;
aircraft.wing.sweep = 50;
aircraft.wing.T2C = .055;   % thickness to chord
aircraft.wing.l = 50;
aircraft.wing.x_c = .24;

aircraft.ht.VolCoeff      = 0.40;
aircraft.ht.AR            = 0.8*aircraft.wing.AR;
aircraft.ht.TaperRatio    = 0.4;
aircraft.ht.leverArm_frac = 0.35;
aircraft.ht.sweep = 30;
aircraft.ht.T2C = .05;
aircraft.ht.x_c = .24;

aircraft.vt.VolCoeff      = 0.04;
aircraft.vt.AR            = 1.8;
aircraft.vt.TaperRatio    = 0.30;
aircraft.vt.leverArm_frac = 0.30;
aircraft.vt.twinTail      = true;
aircraft.vt.sweep = 35;
aircraft.vt.T2C = .05;
aircraft.vt.x_c = .24;

aircraft.ordinance.weight = 4380;

aircraft.avionics.weight = 2500;

% Empty Weight Buildup
aircraft.constants.limitLoad = 8;
aircraft.constants.negLimitLoad = -3;
aircraft.constants.ultLoad = 1.5 * aircraft.constants.limitLoad;
aircraft.constants.maxMach = 1.8;

% CG and Inertia
aircraft.engine.cg.x = 45;
aircraft.cockpit.cg.x = 8;
aircraft.cockpit.weight = 300; % undefined? %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Landing Gear
aircraft.gear.mg.extendedLength = 40;   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% fixxxxxxxxxxxxxxxxxxxx
aircraft.gear.ng.extendedLength = 40;
aircraft.gear.mg.x = 40;
aircraft.gear.ng.x = 8;
aircraft.gear.mg.height = 5;
aircraft.gear.ng.height = 5; 

% Time-Step Mission


% Tolerances
aircraft.weight.tolerance = 15; % GO TO getConfig to uncomment
aircraft.cg.tolerance = 3/12;
aircraft.gear.tolerance = 3/12;

aircraft = dimensionalize_aircraft(aircraft);
aircraft = aeroupdater(aircraft);
% aircraft = EWB(aircraft);

% Ensure required sub-structs exist
if ~isfield(aircraft, 'constants') || ~isstruct(aircraft.constants)
    aircraft.constants = struct();
end

if ~isfield(aircraft, 'weight') || ~isstruct(aircraft.weight)
    aircraft.weight = struct();
end

% Only set defaults if they do not already exist
if ~isfield(aircraft.constants, 'wingLoading')
    aircraft.constants.wingLoading = 115;
    warning("Used hardcoded wingLoading")
end

if ~isfield(aircraft.weight, 'tolerance')
    aircraft.weight.tolerance = 0.06;
    warning("Used hardcoded weight.tolerance")
end


% aircraft = A2A(aircraft);
aircraft = Strike(aircraft);
fprintf('weight: %0.1f lbf', aircraft.weight.total)