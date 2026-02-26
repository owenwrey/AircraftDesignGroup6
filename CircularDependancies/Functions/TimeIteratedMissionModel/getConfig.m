function f = getConfig(aircraft)
%% run "cfg = getConfig();"
% at the beginning of each script, loading all variables in this function 
% into the script as "cgf.exampleVar1", "cfg.exampleVar2", etc.

f.wingLoading = aircraft.constants.wingLoading; % takeoff wing loading [psf]
f.thrust = 58000; % [lb]

f.W.TOguess = 63738; % [lb]
f.W.fuelReq = 16440; % [lb]
f.weightTolerance = aircraft.weight.tolerance;
f.W.crew = 300; % [lb]
f.W.PL.A2A = 2390; % [lb]
f.W.PL.STK = 4380; % [lb]
f.fuelBufferPercent = 0.06;

% Re-evaluate the following using CL/CD and CL^(3/2)/CD optimization
f.cruise.altitude = 30000; % [ft]

f.cruise.distance.out.A2A = 700; % [nmi]
f.cruise.speed.out.A2A = 468; % [KTAS] not yet implemented
f.cruise.distance.in.A2A = 700; % [nmi]
f.cruise.speed.in.A2A = 414; % not yet implemented
f.loiter1.speed.A2A = 212; % KTAS
f.loiter1.time.A2A = 45; % not yet implemented
f.loiter2.time.A2A = 0; % not yet implemented

f.cruise.distance.out.STK = 700; % [nmi]
f.cruise.speed.out.STK = 468; % [KTAS] not yet implemented
f.cruise.distance.in.STK = 700; % [nmi]
f.cruise.speed.in.STK = 414; % not yet implemented
f.loiter1.speed.STK = 212; % KTAS
f.loiter1.time.STK = 45; % not yet implemented
f.loiter2.time.STK = 0; % not yet implemented

f.descentRate = -3000; % [ft/min]

% f.dragPolar.CD0 = 0;
% f.dragPolar.K1 = 0;
% f.dragPolar.K2 = 0;
% f.dragPolar.CDR.A2A = 0;
% f.dragPolar.CDR.STK = 0;
f.CLmax.clean = 1.2;

f.thrustToWeight_idle = 0.05;

f.SFC.idle = 0.5; % [lb/(lb.h)]
f.SFC.takeoff = 1.85;
f.SFC.climb = 0.75;
f.SFC.cruise = 0.7; % estimate
f.SFC.loiter = 0.7;
f.SFC.combat = 1; % reevaluate
f.SFC.weightdrop = 1;
f.SFC.descent = 0.6;
f.SFC.shutdown = 0.5;

end