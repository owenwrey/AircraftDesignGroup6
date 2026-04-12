clc; clear; close all

AR_list = 3:0.5:6;
taper_list = 0.2:0.05:0.5;
sweep_list = 20:2:40;

bestJ = inf;
bestAircraft = [];
bestAR = NaN;
bestTaper = NaN;
bestSweep = NaN;

for AR = AR_list
    for taper = taper_list
        for sweep = sweep_list
            try
                [J, aircraft] = objective_cd0(AR, taper, sweep);

                fprintf('AR=%.2f taper=%.2f sweep=%.2f --> J=%.6f\n', ...
                    AR, taper, sweep, J);

                if J < bestJ
                    bestJ = J;
                    bestAircraft = aircraft;
                    bestAR = AR;
                    bestTaper = taper;
                    bestSweep = sweep;
                end

            catch ME
                fprintf('FAILED: AR=%.2f taper=%.2f sweep=%.2f\n', ...
                    AR, taper, sweep);
                fprintf('Reason: %s\n', ME.message);
            end
        end
    end
end

fprintf('\nBest design found:\n');
fprintf('AR = %.2f\n', bestAR);
fprintf('Taper ratio = %.2f\n', bestTaper);
fprintf('Sweep = %.2f deg\n', bestSweep);
fprintf('Best objective J = %.6f\n', bestJ);