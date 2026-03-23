function Aircraft = landingGear(Aircraft)

%% -------------------------
% Inputs 
% --------------------------
x_cg   = Aircraft.cg.x;                  % [ft] 
x_tail = Aircraft.fuselage.length;       % [ft]
b      = Aircraft.wing.span;             % [ft]
b_half = b/2;                            % [ft]

theta_tipback   = 18;    % [deg]
theta_tailstrike= 12;    % [deg]
theta_min       = 15;    % [deg]
phi_roll        = 5;     % [deg]
clearance_req   = 6/12;  % [ft]
f_nose          = 0.12;  
psi_max         = 63;    % [deg]

z_cg_body = Aircraft.cg.z;          % [ft] CG height above ground

%% -------------------------
% Initial guess
% --------------------------
L = 40/12;        % [ft]
tol = 1e-3;       % [ft]
err = 1;
iter = 0;
iter_max = 100;

%% Governing pitch constraint
theta_design = max([theta_tipback, theta_tailstrike, theta_min]);

%% -------------------------
% Iteration loop
% --------------------------
while err > tol && iter < iter_max
    
    iter = iter + 1;
    L_old = L;
    
    % CG height above ground
    H = z_cg_body;   % already in ft
    
    % Main gear longitudinal location
    a = H * tand(theta_design);   % [ft]
    x_m = x_cg + a;               % [ft]
    
    % Tail distance
    x_t = x_tail - x_m;           % [ft]
    
    % Gear height from tail strike
    L_tail = x_t * tand(theta_tipback);
    
    % Gear height from roll clearance
    L_roll = clearance_req + b_half * sind(phi_roll);
    
    % Controlling gear height
    L = max(L_tail, L_roll);
    
    % Convergence
    err = abs(L - L_old);
end

%% -------------------------
% Nose gear location
% --------------------------
x_n = x_m - (x_m - x_cg)/f_nose;

%% -------------------------
% Main gear track width
% --------------------------
track_half = H / tand(psi_max);
T_main = 2 * track_half;

%% -------------------------
% Store results
% --------------------------
Aircraft.gear.mg.x   = x_m;
Aircraft.gear.ng.x   = x_n;
Aircraft.gear.height   = L;
Aircraft.gear.track    = T_main;


%% Output

fprintf('Iterations                  = %d\n', iter);
fprintf('Landing gear height L       = %.3f ft\n', L);
fprintf('CG height above ground H    = %.3f ft\n', H);
fprintf('Main gear station x_m       = %.3f ft\n', x_m);
fprintf('Nose gear station x_n       = %.3f ft\n', x_n);
fprintf('Tailstrike distance x_t     = %.3f ft\n', x_t);
fprintf('Tailstrike-driven height    = %.3f ft\n', L_tail);
fprintf('Roll-driven height          = %.3f ft\n', L_roll);
fprintf('Main gear half-track        = %.3f ft\n', track_half);
fprintf('Main gear full track        = %.3f ft\n', T_main);

end