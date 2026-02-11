function p = initParams()
%INITPARAMS Constants and unit conversions used by all mission segments.

p.W_S   = 115;          % takeoff wing loading (psf)
p.rho_SL = 0.0023769;   % sea-level density (slug/ft^3)

p.mps2kts = 1.94384;    % m/s -> kt
p.kts2fps = 1/0.59248;  % kt -> ft/s
p.NM2ft   = 6067;       % nmi -> ft
p.ft2m    = 0.305;      % ft -> m

p.Thrust_SLS = 58000;   % lb (reference thrust)

% Thrust lapse model: Altitude in ft -> lapse factor [-]
p.TLapse = griddedInterpolant( ...
    [0, 10000, 20000, 30000, 40000, 50000], ...
    [1, 0.80, 0.60, 0.40, 0.20, 0.15], ...
    'linear','nearest');

% Default SFCs (lb/(lbf*hr))
p.SFC_idle     = 0.50;
p.SFC_takeoff  = 1.85;
p.SFC_climb    = 0.75;
p.SFC_cruise   = 0.70;
p.SFC_descent  = 0.60;
p.SFC_combat   = 1.00;
p.SFC_wtDrop   = 1.00;
p.SFC_loiter   = 0.70;
p.SFC_shutdown = 0.50;

% Segment-specific fixed numbers (edit as needed)
p.cruiseOut_nm = 650;
p.cruiseIn_nm  = 650;

p.cruiseOut_KTAS = 468;
p.cruiseIn_KTAS  = 414;

p.combat_nm   = 100;
p.combat_KTAS = 595;    % ~Mach 0.9 at SL in your original script

p.descent_fpm = -3000;

p.loiter_min  = 20;
p.loiter_KTAS = 212;

% Payload / stores dropped (lb)
p.payloadDrop_lb = 4380;
end
