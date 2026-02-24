%% simpleGearPlacement.m
% Minimal tricycle gear placement calculator + gear height sizing

clc; clear;

%% Inputs (meters, degrees)
x_cg = 6.50;        % CG x-location (m) % this will come from Kyle
z_cg = 1.20;        % CG height above ground (m)n this will come from Kyle

x_tail = 0.80;      % tailstrike point x-location (m)
z_tail = 0.35;      % tailstrike point height above ground (m), CURRENT stance

track_main = 3.20;  % main gear track width (m)

x_m = 7.20;         % chosen main gear x-location (m), must be > x_cg
f_n = 0.10;         % desired nosewheel weight fraction (Rn/W), typical 0.08–0.15

% Requirements
tipback_min_deg  = 15;     % minimum tailstrike margin angle (deg)
overturn_max_deg = 63;     % max overturn angle (deg), 54 if carrier-based

% NEW: Ground attitude + rotation requirement for gear height sizing
theta_req_deg    = 15;     % required rotation/tailstrike angle requirement (deg)
theta_ground_deg = 2;      % desired parked ground pitch angle (deg), nose-up positive

%% 1) Compute nose gear location from static moments
% f_n = (x_m - x_cg)/(x_m - x_n)  ->  x_n = x_m - (x_m - x_cg)/f_n
x_n = x_m - (x_m - x_cg)/f_n;

%% 2) Basic checks (angles from CURRENT stance)
tipback_deg  = atan2d(z_tail, (x_m - x_tail));      % tailstrike margin angle (deg)
overturn_deg = atan2d((track_main/2), z_cg);        % rear-view overturn angle (deg)

wheelbase = x_m - x_n;                              % distance between nose and main gear (m)

%% 3) NEW: Compute vertical gear lengths (heights) from geometry constraints
% We assume changing main gear height lifts the aircraft, increasing tail height by the same amount.
% You can treat this as first-pass sizing.
tail_arm = x_m - x_tail;                            % horizontal distance main gear -> tail point (m)

% Required tail height to meet rotation/tailstrike requirement:
% theta_req = atan(z_tail_req / tail_arm)  -> z_tail_req = tail_arm * tan(theta_req)
z_tail_req = tail_arm * tand(theta_req_deg);        % required tail height above ground (m)

% Additional lift needed compared to current stance:
h_main = max(0, z_tail_req - z_tail);               % main gear height increase needed (m)

% Nose gear height to achieve desired ground pitch angle:
% tan(theta_ground) = (h_nose - h_main)/wheelbase  -> h_nose = h_main + wheelbase*tan(theta_ground)
h_nose = h_main + wheelbase * tand(theta_ground_deg);

% Updated tail height and tailstrike angle after applying h_main
z_tail_new = z_tail + h_main;
tipback_deg_new = atan2d(z_tail_new, tail_arm);

%% 4) Print results + pass/fail
fprintf("\n--- Simple Tricycle Gear Placement ---\n");
fprintf("x_main        = %.3f m\n", x_m);
fprintf("x_nose        = %.3f m\n", x_n);
fprintf("wheelbase     = %.3f m\n", wheelbase);
fprintf("nose frac     = %.3f (%.1f%%)\n", f_n, 100*f_n);

fprintf("\n--- Stability Checks (from CURRENT stance) ---\n");
fprintf("tipback       = %.2f deg  (%s)\n", tipback_deg, passfail(tipback_deg >= tipback_min_deg));
fprintf("overturn      = %.2f deg  (%s)\n", overturn_deg, passfail(overturn_deg <= overturn_max_deg));

fprintf("\n--- Gear Height Sizing (first-pass) ---\n");
fprintf("theta_req     = %.2f deg\n", theta_req_deg);
fprintf("theta_ground  = %.2f deg\n", theta_ground_deg);
fprintf("h_main        = %.3f m  (main gear vertical length/height)\n", h_main);
fprintf("h_nose        = %.3f m  (nose gear vertical length/height)\n", h_nose);
fprintf("tipback NEW   = %.2f deg  (%s)\n", tipback_deg_new, passfail(tipback_deg_new >= tipback_min_deg-1e-6));

%% local helper
function s = passfail(tf)
if tf, s = "PASS"; else, s = "FAIL"; end
end