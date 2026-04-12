function aircraft = runAircraftSizing(AR, taper_ratio, sweep)

addpath(genpath('Functions'));

missionToRun = "Strike";

%% Aircraft struct
aircraft = struct;

% Initial guesses
aircraft.cg.x = 20;
aircraft.cg.y = 0;
aircraft.cg.z = 10;

aircraft.weight.total = 60e3;
aircraft.weight.fuel = 20e3;
aircraft.weight.empty = 40e3;
aircraft.weight.totalOnLanding = 42e3;

%% Constants
aircraft.constants.wingLoading = 102; % [lbf/ft^2]
aircraft.engine.weight = 5000;
aircraft.engine.thrust = 35000;
aircraft.engine.thrustMil = 26000;
aircraft.engine.TSFC = 0.67;

aircraft.constants.fuelVolume = 3500;
aircraft.fuelSys.VI = 0;
aircraft.fuelSys.VP = aircraft.constants.fuelVolume/2;
aircraft.fuelSys.Nt = 4;

%% Fuselage
aircraft.fuselage.length = 48;
aircraft.fuselage.diameter = 6;
aircraft.fuselage.cg.x = 24;
aircraft.fuselage.cg.y = 0;
aircraft.fuselage.cg.z = 10;

%% Wing design variables
aircraft.wing.AR = AR;
aircraft.wing.taper_ratio = taper_ratio;
aircraft.wing.sweep = sweep;

%% Wing fixed variables
aircraft.wing.T2C = 0.055;
aircraft.wing.l = 50;
aircraft.wing.x_c = 0.24;

%% Horizontal tail
aircraft.ht.VolCoeff = 0.40;
aircraft.ht.AR = 4;
aircraft.ht.TaperRatio = 0.4;
aircraft.ht.leverArm_frac = 0.35;
aircraft.ht.sweep = 30;
aircraft.ht.T2C = 0.05;
aircraft.ht.x_c = 0.24;

%% Vertical tail
aircraft.vt.VolCoeff = 0.04;
aircraft.vt.AR = 1.8;
aircraft.vt.TaperRatio = 0.30;
aircraft.vt.leverArm_frac = 0.30;
aircraft.vt.twinTail = true;
aircraft.vt.sweep = 35;
aircraft.vt.T2C = 0.05;
aircraft.vt.x_c = 0.24;

%% Other weights
aircraft.ordinance.weight = 4380;
aircraft.avionics.weight = 2500;
aircraft.engine.cg.x = 45;
aircraft.cockpit.cg.x = 8;
aircraft.cockpit.weight = 300;

%% Loads
aircraft.constants.limitLoad = 8;
aircraft.constants.negLimitLoad = -3;
aircraft.constants.ultLoad = 1.5 * aircraft.constants.limitLoad;
aircraft.constants.maxMach = 1.8;

%% Landing gear
aircraft.gear.mg.extendedLength = 40;
aircraft.gear.ng.extendedLength = 40;
aircraft.gear.mg.x = 40;
aircraft.gear.ng.x = 8;
aircraft.gear.mg.height = 5;
aircraft.gear.ng.height = 5;

%% Tolerances
aircraft.weight.tolerance = 15;
aircraft.cg.tolerance = 3/12;
aircraft.gear.tolerance = 3/12;

%% Iteration loop
exitFlag = false;
iteration = 0;
iterationMax = 1000;

while ~exitFlag && iteration <= iterationMax
    iteration = iteration + 1;
    aircraftOld = aircraft;

    aircraft = dimensionalize_aircraft(aircraft);
    aircraft = aeroupdater(aircraft);
    aircraft = EWB(aircraft);
    aircraft = CgInertiaCalc(aircraft);
    aircraft = roskam_inertia(aircraft);
    aircraft = landingGear(aircraft);

    if abs(aircraftOld.gear.mg.x - aircraft.gear.mg.x) > aircraft.gear.tolerance && ...
       abs(aircraftOld.gear.ng.x - aircraft.gear.ng.x) > aircraft.gear.tolerance
        continue;
    end

    if abs(aircraftOld.weight.total - aircraft.weight.total) > aircraft.weight.tolerance
        continue;
    end

    aircraft = TIMESTEP_CONVERGENCE_MASTER(aircraft, missionToRun);

    if abs(aircraftOld.weight.total - aircraft.weight.total) < aircraft.weight.tolerance && ...
       abs(aircraftOld.cg.x - aircraft.cg.x) < aircraft.cg.tolerance
        exitFlag = true;
    end

    aircraft.constants.fuelVolume = aircraft.weight.fuel/6.7;
    aircraft.fuelSys.VP = aircraft.constants.fuelVolume/2;
    aircraft.constants.thrustToWeight_TO.AB = 2*aircraft.engine.thrust/aircraft.weight.total;
    aircraft.constants.thrustToWeight_TO.mil = 2*aircraft.engine.thrustMil/aircraft.weight.total;
    aircraft.constants.EWF = aircraft.weight.empty/aircraft.weight.total;
end

if iteration > iterationMax
    error('Aircraft sizing did not converge.');
end