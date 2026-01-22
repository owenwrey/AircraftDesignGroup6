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
f.cruise.altitude = 25000; % [ft]
f.cruise.distance.out.A2A = 0;
f.cruise.speed.out.A2A = 0; % [KTAS] not yet implemented
f.cruise.distance.in.A2A = 850; % [nmi]
f.cruise.speed.in.A2A = 635; % not yet implemented
f.loiter1.speed.A2A = 481; % KTAS
f.loiter1.time.A2A = 45; % not yet implemented
f.loiter2.time.A2A = 0; % not yet implemented

% f.dragPolar.CD0 = 0;
% f.dragPolar.K1 = 0;
% f.dragPolar.K2 = 0;
% f.dragPolar.CDR = 0;

f.wingLoading = 130; % takeoff wing loading [psf]
f.thrust = 44000; % [lb]

f.thrustToWeight_idle = 0.05;

f.SFC.idle = 0.6; % [lb/(lb.h)]
f.SFC.takeoff = 1.85;
f.SFC.climb = 0.75;
f.SFC.cruise = 0.85; % estimate
f.SFC.loiter = 0.75;
f.SFC.combat = 0.85; % reevaluate
f.SFC.weightdrop = 0.85;
f.SFC.descent = 0.60;
f.SFC.shutdown = 0.6;

end