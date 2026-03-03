%% TurnPerformance_PsLines_20kft.m
% Makes a turn performance chart like TurnPerformanceChatGPT.m, but:
%   - fixed altitude = 20,000 ft
%   - each curve corresponds to a target excess power / climb rate (Ps = hdot)
%
% Assumptions:
%   - coordinated, level turn kinematics for psi-dot and turn radius
%   - Ps(V,n) = (T - D(V,n))*V/W  [ft/s]  (same as climb script style)
%   - Drag uses getCD(CL,M) from TurnPerformanceChatGPT.m

clear variables; close all; clc;

%% Unit conversions / constants
fps2kn  = 0.592484;
nmi2ft  = 6076.12;
gc      = 32.174;      % [lbm/slug] (used as ft/s^2 conversion in your scripts)

%% User settings
alt_ft = 20e3;         % fixed altitude [ft]
numV   = 160;          % number of speed points along each curve

% Pick lines by climb rate (equivalently Ps). Edit to taste.
climbRate_fpm_list = [- 500 0 500 1000 1500 2000 3000 4000];   % [ft/min]
PsTarget_fps_list  = climbRate_fpm_list / 60;            % [ft/s] (Ps = hdot)

% Background contour settings (same style as TurnPerformanceChatGPT.m)
R_levels_nmi = [0.2 0.3 0.5 0.8 1 1.5 2 3 5 8 12];
n_levels_g   = [2 3 4 5 6 7 8];

%% Aircraft parameters (match your existing scripts)
W0   = 84127;          % [lbf]
Sref = W0/115;         % [ft^2]
W    = 72614;          % [lbf] (you used this in the other scripts)
T_SL = 58000;          % [lbf] max thrust sea level static

CLmax       = 1.2;
nMaxStruct  = 8.0;     % [g]

%% Atmosphere at 20k ft
[~, a_fps, ~, rho_slugft3] = atmosEnglish(alt_ft, 'ft', 'slug', 'R');

%% Thrust at 20k ft (same lapse model as TurnPerformanceChatGPT.m)
TLapse = griddedInterpolant( ...
    [0, 10000, 20000, 30000, 40000, 50000], ...
    [1, 0.80, 0.60, 0.40, 0.20, 0.15], ...
    'linear','nearest');

alpha_thrust = TLapse(alt_ft);
T_available  = T_SL * alpha_thrust;   % [lbf] (constant vs Mach here, consistent w/ turn script)

%% Build speed grid (use 1g stall to max-thrust-limited Vmax like TurnPerformanceChatGPT.m)
Vs1g_fps = sqrt( (2/rho_slugft3) * (W/Sref) * (1/CLmax) );

% Find Vmax where excess thrust at n=1 goes to zero: T - D(V, n=1) = 0
Pexc_func = @(V) T_available - 0.5*rho_slugft3*V.^2*Sref .* getCD( ...
    2*W./(rho_slugft3*V.^2*Sref), ...  % CL at n=1
    V./a_fps );                        % Mach

options = optimset('TolX', 1e-3);
Vmax_fps = fzero(Pexc_func, [Vs1g_fps, 6000], options);

V_fps  = linspace(Vs1g_fps, Vmax_fps, numV);
V_KTAS = V_fps * fps2kn;

q_psf  = 0.5*rho_slugft3 .* V_fps.^2;

%% Precompute instantaneous (stall/struct) n_max at each speed
nMaxInst = q_psf * CLmax / (W/Sref);
nMaxInst = min(nMaxInst, nMaxStruct);
nMaxInst(1) = 1.0; % at stall end of the grid, enforce 1g

%% For each Ps target, compute max-n satisfying Ps(V,n) >= PsTarget
numPs = numel(PsTarget_fps_list);
nMax_forPs      = NaN(numPs, numV);
psiDot_forPs    = NaN(numPs, numV);

for iPs = 1:numPs
    PsTarget = PsTarget_fps_list(iPs);

    for j = 1:numV
        V = V_fps(j);
        q = q_psf(j);
        M = V/a_fps;

        % Helper: Ps(V,n)
        Ps_of_n = @(n) (T_available - q*Sref*getCD( n*W/(q*Sref), M )) * V / W;

        % Feasibility checks
        Ps_at_1g = Ps_of_n(1.0);
        if Ps_at_1g < PsTarget
            % Can't even meet Ps at level (n=1) at this V
            continue;
        end

        n_hi = nMaxInst(j);
        Ps_at_hi = Ps_of_n(n_hi);

        if Ps_at_hi >= PsTarget
            % Even the instantaneous limit can meet PsTarget; use it
            n_sol = n_hi;
        else
            % Solve Ps(n) = PsTarget via bisection on n in [1, n_hi]
            n_lo = 1.0;

            for k = 1:50
                n_mid = 0.5*(n_lo + n_hi);
                Ps_mid = Ps_of_n(n_mid);

                if Ps_mid >= PsTarget
                    n_lo = n_mid;   % can still meet Ps, try higher n
                else
                    n_hi = n_mid;   % too much drag, reduce n
                end
            end
            n_sol = n_lo; % max n that still meets PsTarget
        end

        nMax_forPs(iPs,j)   = n_sol;
        psiDot_forPs(iPs,j) = gc*sqrt(max(n_sol^2 - 1, 0))/V * (180/pi);
    end
end

%% Build background contour grid in (V, psiDot): iso-n and iso-R
Vg_kn = linspace(min(V_KTAS), max(V_KTAS), 450);
ng    = linspace(1.0, nMaxStruct, 450);

[VG_kn, NG] = meshgrid(Vg_kn, ng);

VG_fps = VG_kn / fps2kn;

psiDotG = (gc .* sqrt(NG.^2 - 1) ./ VG_fps) * (180/pi);       % [deg/s]
RG_nmi  = (VG_fps.^2 ./ (gc .* sqrt(NG.^2 - 1))) / nmi2ft;    % [nmi]

psiDotG(NG <= 1.0001) = NaN;
RG_nmi (NG <= 1.0001) = NaN;

%% Plot
figure; hold on; grid on; box on;

% iso-R thin
[CR, hR] = contour(VG_kn, psiDotG, RG_nmi, R_levels_nmi, 'LineWidth', 0.9);
clabel(CR, hR, 'FontSize', 9);

% iso-n thicker
[Cn, hn] = contour(VG_kn, psiDotG, NG, n_levels_g, 'LineWidth', 0.9);
clabel(Cn, hn, 'FontSize', 9);

% Overlay Ps lines at 20k ft
lw = 1.8;
for iPs = 1:numPs
    plot(V_KTAS, psiDot_forPs(iPs,:), '-', 'LineWidth', lw);
end

xlabel('KTAS', 'Interpreter','latex');
ylabel('Turn Rate, $\dot{\Psi}$ (deg/s)', 'Interpreter','latex');
title(sprintf('Turn Performance at %.0f ft: Lines of Excess Power / Climb Rate', alt_ft), 'Interpreter','latex');

% y-limit similar to your existing script (keep readable)
yMax = max(psiDot_forPs(:), [], 'omitnan');
ylim([0, min(1.05*yMax, 40)]);

% Legend text for Ps lines
psLabels = strings(1,numPs);
for iPs = 1:numPs
    psLabels(iPs) = sprintf('Climb rate = %.0f fpm', climbRate_fpm_list(iPs));
end

legend([ ...
    "Turn Radius (nmi)", "Load Factor (g)", psLabels], ...
    'Location','northeastoutside', 'Interpreter','latex');

%% ---------------- FUNCTIONS ----------------
function [T_Eng, SOS_fps, p_psf, rhoEng] = atmosEnglish(h, altUnits, densityUnits, tempUnits)
% atmosisa conversion to English Engineering units
% input : (h, altUnits, densityUnits, tempUnits)
% output: [T_Eng, SOS_fps, p_psf, rhoEng]
% defaults units are feet (input), slugs/ft^3, and Rankine

ft2m = 0.3048;
K2Rankine = 9/5;
mps2fps = 1/ft2m;
Pa2psf = 0.02088543427;
kg_m3_2lbm_cuft = 0.06242796;

gc = 32.174; % [lbm/slug]

if nargin < 2 || isempty(altUnits), altUnits = 'ft'; end
if strcmpi(altUnits,'m') || strcmpi(altUnits,'meter') || strcmpi(altUnits,'meters')
    h_m = h;
elseif strcmpi(altUnits,'ft') || strcmpi(altUnits,'feet') || strcmpi(altUnits,'foot')
    h_m = h * ft2m;
else
    error('Altitude units not recognized.')
end

[T_K, SOS_mps, p_Pa, rho_kg_m3, ~, ~] = atmosisa(h_m);

T_Ra     = K2Rankine*T_K;
SOS_fps  = mps2fps*SOS_mps;
p_psf    = Pa2psf*p_Pa;
rho_lbm  = kg_m3_2lbm_cuft*rho_kg_m3;
rho_slug = rho_lbm / gc;

if nargin < 3 || isempty(densityUnits), densityUnits = 'slug'; end
if strcmpi(densityUnits,'lbm') || strcmpi(densityUnits,'lb') || strcmpi(densityUnits,'pounds') || strcmpi(densityUnits,'lbs')
    rhoEng = rho_lbm;
elseif strcmpi(densityUnits,'slug') || strcmpi(densityUnits,'slugs') || strcmpi(densityUnits,'sl') || strcmpi(densityUnits,'slug/ft^3')
    rhoEng = rho_slug;
else
    error('Density units not recognized.')
end

if nargin < 4 || isempty(tempUnits), tempUnits = 'Rankine'; end
if strcmpi(tempUnits,'Ra') || strcmpi(tempUnits,'R') || strcmpi(tempUnits,'Rankine')
    T_Eng = T_Ra;
elseif strcmpi(tempUnits,'K') || strcmpi(tempUnits,'Kelvin')
    T_Eng = T_K;
else
    error('Temperature units not recognized.')
end
end

function CD = getCD(CL, M)
% Drag polar (copied from TurnPerformanceChatGPT.m)
CD0_sub = 0.02;
K2_sub  = 0.025;

CD0_super = 0.022;
K2_super  = 0.025;

CDR = zeros(size(M));

machCurve = [0        0.85    0.95                               1.05       2];
CD0curve  = [CD0_sub  CD0_sub 1.125*(CD0_super-CD0_sub)+CD0_sub  CD0_super  0.95*CD0_super];
K2curve   = [K2_sub   K2_sub  1.125*(K2_super-K2_sub)+K2_sub     K2_super   K2_super];

CD0 = interp1(machCurve, CD0curve, M, "linear", "extrap");
K2  = interp1(machCurve, K2curve,  M, "linear", "extrap");

CD = CD0 + K2 .* (CL.^2) + CDR;
end