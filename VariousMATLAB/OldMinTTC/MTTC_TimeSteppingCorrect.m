% Minimum Time to Climb
clear; clc; close all;

% sea level to 46,250 ft

%% Update checklist, MAKE SURE TO DO THIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% things to update include: MTOW, weight fraction, wing loading, max alt,
% thrust, drag values, TSFC, min/max speed in speed array, EnHt start and end, 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Inputs
TLapse = griddedInterpolant([0, 10000, 20000, 30000, 40000, 50000], ...
                            [1, 0.80, 0.60, 0.40, 0.20, 0.15], ...
                            'linear','nearest');

rho_SL = 0.0023769;
mps2kts = 1.94384;      % meters per sec to knots
kts2fps = 1/0.59248;    % knots to feet per sec
NM2ft = 6067;           % nautical miles to feet
ft2m = 0.305;           % feet to meters

W0 = 50000;             % MTOW (lbm)
WtFrac = 1;             % weight fraction
W_S = 112;              % takeoff wing loading (psf)
S = W0/W_S;             % wing area ft^2
maxAlt = 49000;         % service ceiling, ft
fullThrust = 58000;     % lb

CD0 = 0.02;            
K1 = 0;      
K2 = 0.05;   
CDR = 0;     

istart = 1;
iend = 250;

SFC_climb = 0.75;
THROT = 1;

Alt = linspace(0,maxAlt, iend)';   % altitude array
[T, a, P, rho] = atmosisa(Alt*ft2m);
rho = rho*0.00194032; % to slug/ft^3
thrustLapse = TLapse(Alt);

KTAS = linspace(160,1000,iend)';                % generate speed array

EnHt_start = (KTAS(1)*kts2fps).^2/(2*32.2);
EnHt_end   = maxAlt + (600*kts2fps).^2/(2*32.2);

EnHt = linspace(EnHt_start,EnHt_end,(iend-2))';

%% Build the grids

[ALT_grid, KTAS_grid] = meshgrid(Alt, KTAS);    % generate speed and alt grid
V_grid = KTAS_grid * kts2fps;

q_grid = 0.5 .* (rho') .* (V_grid.^2);          % generate q grids

CL_grid = W0 ./ (q_grid .* S);                  % generate grids for CL and CD
CD_grid = CD0 + K1 * CL_grid + K2 * CL_grid.^2;
% code here assumes weight is constant throughout climb. Later model should
% have time stepping and calculate fuel burn for more accurate modeling

Drag_grid = q_grid .* S .* CD_grid;
Thrust_grid = fullThrust .* thrustLapse';

Ps_grid = (Thrust_grid - Drag_grid) .* V_grid ./ W0;

% Energy height map
En_grid = ALT_grid + (V_grid.^2)/(2*32.2);

%% Best Excess Power Logic

tol = 2500;   % wide tolerance for matching energy height

bestPs  = zeros(length(EnHt),1);
bestAlt = zeros(length(EnHt),1);
bestV   = zeros(length(EnHt),1);

Vmax_kts = 900;                 % max allowed KTAS
Vmax = Vmax_kts * kts2fps;      % ft/s

for i = 1:length(EnHt)
    mask = abs(En_grid - EnHt(i)) < tol & V_grid <= Vmax;
    Ps_candidates = Ps_grid(mask);

    if isempty(Ps_candidates)
        bestPs(i)  = NaN;
        bestAlt(i) = NaN;
        bestV(i)   = NaN;
        continue
    end

    [bestPs(i), idx] = max(Ps_candidates);

    ALT_candidates = ALT_grid(mask);
    V_candidates   = V_grid(mask);

    bestAlt(i) = ALT_candidates(idx);
    bestV(i)   = V_candidates(idx);
end

%% Remove NaNs from the climb path

good = ~isnan(bestPs);
bestPs  = bestPs(good);
bestAlt = bestAlt(good);
bestV   = bestV(good);
EnHt    = EnHt(good);

%% Force Initial State: 0 ft, 200 TKAS

initAlt = 0;
initV_kts = 200;
initV = initV_kts*kts2fps;

q_init = 0.5*rho(1)*initV^2;

CL_init = W0/(q_init*S);
CD_init = CD0 + K1*CL_init + K2*CL_init^2;

D_init = q_init*S*CD_init;
T_init = fullThrust*thrustLapse(1);

Ps_init = (T_init - D_init)*initV/W0;

EnHt_init = initAlt + initV^2/(2*32.2);

bestAlt = [initAlt; bestAlt];
bestV   = [initV;   bestV];
bestPs  = [Ps_init; bestPs];
EnHt    = [EnHt_init; EnHt];


%% Force Final State: 46250 ft, 600 KTAS
finalAlt    = maxAlt;     % 46,250 ft
finalV_kts  = 600;        % 600 KTAS
finalV      = finalV_kts * kts2fps;

q_final = 0.5 * rho(end) * finalV^2;
CL_final = W0 / (q_final * S);
CD_final = CD0 + K2 * CL_final^2;

D_final = q_final * S * CD_final;
T_final = fullThrust * thrustLapse(end);
Ps_final = (T_final - D_final) * finalV / W0;

EnHt_final = finalAlt + finalV^2/(2*32.2);

% Append the mandatory endpoint
bestAlt(end+1) = finalAlt;
bestV(end+1)   = finalV;
bestPs(end+1)  = Ps_final;
EnHt(end+1)    = EnHt_final;

% Sort (important)
[EnHt, idx] = sort(EnHt);
bestAlt = bestAlt(idx);
bestV   = bestV(idx);
bestPs  = bestPs(idx);

% ---- REMOVE DUPLICATE ENERGY HEIGHTS ----
[EnHt, uniqueIdx] = unique(EnHt, 'stable');

bestAlt = bestAlt(uniqueIdx);
bestV   = bestV(uniqueIdx);
bestPs  = bestPs(uniqueIdx);

%% TIME MARCH SIMULATION (FIXED dt)

g = 32.2;                 % ft/s^2
dt = 1.0;                 % time step (sec) — small for stability

% Initial state
t = 0;
Alt = initAlt;
V   = initV;
W   = W0;
He  = Alt + V^2/(2*g);

% Storage arrays
t_hist   = t;
Alt_hist = Alt;
V_hist   = V;
W_hist   = W;

% Convert TSFC to 1/sec (from 1/hr)
TSFC = SFC_climb / 3600;

while Alt < finalAlt

    % Interpolate optimal Ps at current energy height
    Ps_cmd = interp1(EnHt, bestPs, He, 'linear', 'extrap');

    % Atmospheric properties at current altitude
    [~,~,~,rho_now] = atmosisa(Alt*ft2m);
    rho_now = rho_now * 0.00194032;

    thrust_now = fullThrust * TLapse(Alt);
    
    q = 0.5 * rho_now * V^2;
    CL = W/(q*S);
    CD = CD0 + K1*CL + K2*CL^2;
    D = q*S*CD;

    Ps_actual = (thrust_now - D)*V/W;

    % Use commanded Ps but limit to actual available
    Ps = min(Ps_cmd, Ps_actual);

    % Energy integration
    He_dot = Ps;
    He_next = He + He_dot*dt;

    if He_next >= EnHt_final
        % compute fractional time step to land exactly on target
        dt_final = (EnHt_final - He) / He_dot;
        
        He = EnHt_final;
        t  = t + dt_final;
        
        % update V and Alt exactly at final energy
        V = interp1(EnHt, bestV, He, 'linear', 'extrap');
        Alt = He - V^2/(2*g);
        
        % burn fuel during fractional step
        mdot = TSFC * thrust_now;
        W = W - mdot*dt_final;
        
        % store final state
        t_hist(end+1,1)   = t;
        Alt_hist(end+1,1) = Alt;
        V_hist(end+1,1)   = V;
        W_hist(end+1,1)   = W;
        
        break
    else
        He = He_next;
    end


    % Assume optimal schedule splits energy correctly
    % Recompute V and Alt from He and commanded V schedule
    V_cmd = interp1(EnHt, bestV, He, 'linear', 'extrap');
    V_cmd = max(min(V_cmd, max(bestV)), min(bestV));

    V = V_cmd;
    Alt = He - V^2/(2*g);

    % Prevent negative altitude
    Alt = max(0, Alt);

    % Fuel burn
    mdot = TSFC * thrust_now;     % lb/sec
    W = W - mdot*dt;
    W = max(W, 0.7*W0);           % safety floor (avoid negative weight)

    % Time update
    t = t + dt;

    % Store history
    t_hist(end+1,1)   = t;
    Alt_hist(end+1,1) = Alt;
    V_hist(end+1,1)   = V;
    W_hist(end+1,1)   = W;

    
end

while V < finalV
    
    [~,~,~,rho_now] = atmosisa(Alt*ft2m);
    rho_now = rho_now * 0.00194032;

    thrust_now = fullThrust * TLapse(Alt);
    
    q = 0.5 * rho_now * V^2;
    CL = W/(q*S);
    CD = CD0 + K1*CL + K2*CL^2;
    D = q*S*CD;

    accel = (thrust_now - D)/W * g;   % ft/s^2
    
    V = V + accel*dt;

    mdot = TSFC * thrust_now;
    W = W - mdot*dt;

    t = t + dt;

    t_hist(end+1,1)   = t;
    Alt_hist(end+1,1) = Alt;
    V_hist(end+1,1)   = V;
    W_hist(end+1,1)   = W;
end

TTC = t/60;

%% ================================
%% Plotting
%% ================================

figure(1)
plot(V_hist/kts2fps, Alt_hist,'b','LineWidth',2)
xlabel('KTAS')
ylabel('Altitude (ft)')
title('Minimum-Time Climb Schedule (Time Marched)')
grid on

figure(2)
plot(t_hist, Alt_hist,'b','LineWidth',2)
xlabel('Time (s)')
ylabel('Altitude (ft)')
title('Climb Profile')
grid on

figure(3)
plot(t_hist, W_hist,'LineWidth',2)
xlabel('Time (s)')
ylabel('Weight (lb)')
title('Fuel Burn During Climb')
grid on

fprintf('Total Time to Climb = %.2f minutes\n', TTC);
fprintf('Fuel Burned = %.0f lb\n', W0 - W);