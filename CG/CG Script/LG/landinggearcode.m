%% LANDING GEAR LOCATION SOLVER
% Determines main gear & nose gear placement based on:
% - CG location
% - MAC
% - tipback angle requirement (>15 deg)
% - tipover angle requirement (>25 deg)
% - desired nose gear load fraction (10–20%)
% - rotation angle requirement (~10–15 deg)

clear; clc;

% Inputs (edit these for YOUR jet)
CGx = 23.0;      % ft (CG x-location from nose)
CGz = 4.0;       % ft (CG height above ground)
MAC  = 6.128;    % ft
TailHeight = 14; % ft (vertical distance tail sits above ground)
Track = 12;      % ft (distance between main gear legs)
NoseGearLoadFrac = 0.15; % desired W_nose / W_total

% Sweep possible main-gear locations aft of CG 
MG_range = CGx + 0.10*MAC : 0.01 : CGx + 0.25*MAC;

valid = [];

for MGx = MG_range

    % 1) Tip-back angle check
    horiz = MGx - CGx;              % horizontal CG-to-MG distance
    theta_tipback = atan(CGz / horiz) * 180/pi;

    if theta_tipback < 15
        continue; % FAIL → not stable on ground
    end

    % 2) Tip-over angle check
    theta_tipover = atan((Track/2) / CGz) * 180/pi;

    if theta_tipover < 25
        continue; % FAIL → plane tips sideways
    end

    % 3) Compute nose gear location for desired load fraction
    % Solve for NGx using static moment equilibrium:
    % Wn = W * (MGx - CGx) / (MGx - NGx)
    % → NGx = MGx - (MGx - CGx)/NoseGearLoadFrac

    NGx = MGx - (MGx - CGx) / NoseGearLoadFrac;

    if NGx <= 0
        continue; % invalid, nose gear cannot be behind the nose
    end

    % 4) Rotation angle check
    Wheelbase = MGx - NGx;
    theta_rotate = atan(TailHeight / Wheelbase) * 180/pi;

    if theta_rotate < 10 || theta_rotate > 15
        continue; % rotation angle outside ideal range
    end

    % If all conditions passed → store solution
    valid = [valid; MGx, NGx, theta_tipback, theta_tipover, theta_rotate];

end

% Display results
if isempty(valid)
    disp("No valid landing-gear geometry found. Loosen constraints.");
else
    T = array2table(valid, ...
        'VariableNames', {'MainGearX', 'NoseGearX', 'TipbackDeg', 'TipoverDeg', 'RotationDeg'});
    disp(T)
end
