function f = getConfig()
%% run "cfg = getConfig();"
% at the beginning of each script, loading all variables in this function 
% into the script as "cgf.exampleVar1", "cfg.exampleVar2", etc.

f.W.TOguess = 63738;
f.W.fuelReq = 16440;
f.weightTolerance = 2; % [lb]
f.W.crew = 300;
f.W.PL.A2A = 2390;
f.fuelBufferPercent = 0.06;


% Re-evaluate the following using CL/CD and CL^(3/2)/CD optimization
f.cruise.altitude = 30000; % [ft]
f.cruise.distance.out.A2A = 700;
f.cruise.speed.out.A2A = 468; % [KTAS] not yet implemented
f.cruise.distance.in.A2A = 700; % [nmi]
f.cruise.speed.in.A2A = 414; % not yet implemented
f.loiter1.speed.A2A = 212; % KTAS
f.loiter1.time.A2A = 45; % not yet implemented
f.loiter2.time.A2A = 0; % not yet implemented

% f.dragPolar.CD0 = 0;
% f.dragPolar.K1 = 0;
% f.dragPolar.K2 = 0;
% f.dragPolar.CDR = 0;

f.wingLoading = 115; % takeoff wing loading [psf]
f.thrust = 58000; % [lb]

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