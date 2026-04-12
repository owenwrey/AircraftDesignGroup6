function [J, aircraft] = objective_cd0(AR, taper_ratio, sweep)

aircraft = runAircraftSizing(AR, taper_ratio, sweep);

cd0 = aircraft.aero.cd0_strike(30000, 0.9);
W   = aircraft.weight.total;

Wmax = 65000;   % set this to your acceptable upper bound

if W > Wmax
    J = 1e6 + 100*(W - Wmax);
else
    J = cd0;
end