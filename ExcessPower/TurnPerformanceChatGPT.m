%% MASTER run file
% NOT FULLY GENERALIZED
clear variables; close all; clc;

%% Unit conversions
fps2kn = 0.592484;
nmi2ft = 6076.12; % naut. miles to feet
gc = 32.174; % [lbm/slug]

%% Script Settings
% Number of velos in array
numV = 100;
% Get atmos conditions
alt = [0, 15e3, 30e3]; % [ft]
alt = 0:10e3:40e3;
% alt = 20e3;
[~, a, ~, rho] = atmosEnglish(alt, 'ft', 'slug');
numAlt = length(alt);

%% Aircraft Characteristics
W0 = 84127; % Aircraft Weight [lb_f]
Sref = W0/115; % Wing Area [ft^2]
W = 72614;
% W = W0;
T_SL = 58000; % Max Thrust SLS [lb_f]

CLmax = 1.2;

nMaxStruct = 8.0; % [g]

%% Thrust array
TLapse = griddedInterpolant([0, 10000, 20000, 30000, 40000, 50000], [1, 0.80, 0.60, ...
0.40, 0.20, 0.15],'linear','nearest');
alpha_thrust = TLapse(alt); % thrust lapse at 0, 15k, and 30k ft
T = T_SL*alpha_thrust; % row vector of thrusts

% 1 g stall speed
Vs1g = sqrt( (2./(rho)).* (W/Sref)*(1/CLmax) );
% Set Vmin
Vmin = Vs1g;

%% Initialize matrices
V = zeros(numAlt, numV);
q = zeros(numAlt, numV);

nMaxInst  = zeros(numAlt, numV);
nMaxSus = zeros(numAlt, numV);

phiInst  = zeros(numAlt, numV);
phiSus   = zeros(numAlt, numV);

psiDotInst = zeros(numAlt, numV);
psiDotSus  = zeros(numAlt, numV);

R          = zeros(numAlt, numV);

CL         = zeros(numAlt, numV);

%% Find Vmax

for i = 1:numAlt
    
    Vlow = Vs1g(i);
    Vhigh = 6000;
    
    % Excess thrust function (positive when thrust > drag)
    Pexc_func = @(V_var) T(i) - 0.5 .* rho(i) .* V_var.^2 .* Sref .* getCD(2.*W./(rho(i).*V_var.^2.*Sref), V_var./a(i) );
    
    % Use fzero to find where excess thrust = 0
    options = optimset('TolX', 1e-3);
    try
        Vmax(i) = fzero(Pexc_func, [Vlow, Vhigh], options);
    catch ME
        warning('fzero failed for altitude index %d: %s', i, ME.message);
        % Fallback: use upper bound or NaN
        Vmax(i) = NaN;
    end
    
    V(i,:) = linspace(Vs1g(i), Vmax(i), numV); % ignore error
    q(i,:) = 0.5*rho(i)*V(i,:).^2; % get q matrix

    %% Compute Loads
    
    % max load factor before stall

    for j1 = 1:numV
        % Calculate instantaneous n_max
        nMaxInst(i,j1) = q(i,j1).*CLmax./(W/Sref);

        % make sure its not over structural
        if nMaxInst(i,j1) > nMaxStruct
            nMaxInst(i,j1) = nMaxStruct;
        end

        nLow = 1;
        nHigh = nMaxInst(i,j1);
        if nHigh <= nLow
            nMaxSus(i,j1) = NaN;
            continue;
        end
        
        for k = 1:50
            nMid = (nHigh+nLow)/2;
            CL(i,j1) = nMid*W/(q(i,j1)*Sref);

            if CL(i,j1) > CLmax
                nHigh = nMid;
                continue;
            end
            CD = getCD(CL(i,j1), V(j1)/a(i));
            D = q(i,j1)*Sref*CD;
            % D = q(i,j)*Sref*(CD0 + K1*CL(i,j) + K2*CL(i,j)^2);
            
            if T(i) - D > 0
                nLow = nMid;
            else
                nHigh = nMid;
            end
        end % end for k

        nMaxSus(i,j1) = nMid;
        nMaxSus(:,1) = 1.00; % nmax at stall must be 1
    end % end for j1

    %% Compute phi, R, psiDot

    nMaxSus(i,numV) = floor(nMaxSus(i,numV));
    for j2 = 1:numV
        phiSus(i,j2) = acos(1/nMaxSus(i,j2))*180/pi;
        phiInst(i,j2) = real(acos(1/nMaxInst(i,j2))*180/pi); % small imaginary component for some reason

        R(i,j2) = V(i,j2)^2/(gc*sqrt(nMaxSus(i,j2)^2-1))/nmi2ft;

        psiDotSus(i,j2) =  gc*sqrt(nMaxSus(i,j2)^2-1)/V(i,j2)*180/pi;
        psiDotInst(i,j2) = real(gc*sqrt(nMaxInst(i,j2)^2-1)/V(i,j2))*180/pi;
    end
end

%% Plot Combined Turn Performance Chart (V vs PsiDot with iso-n and iso-R)

% convert to knots
V_KTAS = fps2kn * V;   % [KTAS]

figure; hold on; grid on; box on;

% ---------------------------
% 1) Build a background grid for contouring in (V, psidot)
% ---------------------------
% Choose a single V range that spans all altitudes
Vmin_all = min(V_KTAS(:));
Vmax_all = max(V_KTAS(:));

Vg = linspace(Vmin_all, Vmax_all, 450);          % knots
ng = linspace(1.0, nMaxStruct, 450);             % g's (>=1)

[VG, NG] = meshgrid(Vg, ng);                     % grid in (V, n)

% Convert knots -> ft/s for the turn relations (your gc is in ft/s^2)
VG_fps = VG /fps2kn;

% Level coordinated turn relations
psiDotG = (gc .* sqrt(NG.^2 - 1) ./ VG_fps) * (180/pi);    % deg/s
RG_nmi  = (VG_fps.^2 ./ (gc .* sqrt(NG.^2 - 1))) / nmi2ft; % nmi

% Mask the n=1 singular line and any crazy-high rates at low V
psiDotG(NG <= 1.0001) = NaN;
RG_nmi (NG <= 1.0001) = NaN;

% ---------------------------
% 2) Contours: iso-R (thin) and iso-n (thicker)
% ---------------------------
R_levels = [0.2 0.3 0.5 0.8 1 1.5 2 3 5 8 12];   % nmi (tune to taste)
n_levels = [2 3 4 5 6 7 8];                      % g's

% iso-R: thin lines
[CR, hR] = contour(VG, psiDotG, RG_nmi, R_levels, 'LineWidth', 0.9, 'Color', [0.7 0.7 0.9]);
clabel(CR, hR, 'FontSize', 9);

% iso-n: thicker lines
[Cn, hn] = contour(VG, psiDotG, NG, n_levels, 'LineWidth', 0.9, 'Color', [0.7 0.5 0.5]);
clabel(Cn, hn, 'FontSize', 9);

% ---------------------------
% 3) Overlay your sustained / instantaneous boundaries
% ---------------------------
colors = jet(numAlt);

lw = 1.7;

for i = 1:numAlt
    % Sustained boundary
    plot(V_KTAS(i,:), psiDotSus(i,:), '-', ...
        'Color', colors(i,:), 'LineWidth', lw);

    % Instantaneous boundary
    plot(V_KTAS(i,:), psiDotInst(i,:), '--', ...
        'Color', colors(i,:), 'LineWidth', lw);
end

% ---------------------------
% 4) Labels / limits / legend
% ---------------------------
xlabel('KTAS', Interpreter='latex');
ylabel('Turn Rate, $\dot{\Psi}$ (deg/s)', Interpreter='latex');
title('Turn Performance', Interpreter='latex');

% Make y-limits sane (optional but usually necessary)
% Use max of instantaneous data but clamp to something readable
yMax = max(psiDotInst(:), [], 'omitnan');
ylim([0, min(yMax*1.05, 40)]);  % change 40 if you want more/less

legend( ...
    'Turn Radius (nmi)', 'Load Factor (g)', ...
    '0 ft Sustained','0 ft Instantaneous', ...
    '10k ft Sustained','10k ft Instantaneous', ...
    '20k ft Sustained','20k ft Instantaneous', ...
    '30k ft Sustained','30k ft Instantaneous', ...
    '40k ft Sustained','40k ft Instantaneous', ...
    '50k ft Sustained','50k ft Instantaneous', ...
    'Location','northeastoutside', Interpreter='latex');


%% Functions 

function [T_Eng, SOS_fps, p_psf, rhoEng] = atmosEnglish(h, altUnits, densityUnits, tempUnits)
% atmosisa conversion to English Engineering units
% input : (h, altUnits, densityUnits, tempUnits)
% output: [T_Eng, SOS_fps, p_psf, rhoEng]
% defaults units are feet (input), slugs/ft^3, and Rankine

% unit conversions
ft2m = 0.3048;
K2Rankine = 9/5;
mps2fps = 1/ft2m;
Pa2psf = 0.02088543427;
kg_m3_2lbm_cuft = 0.06242796;

gc = 32.174; % [lbm/slug]

% --- altitude units check ---
if nargin < 2 || isempty(altUnits)
    altUnits = 'ft'; % default to feet if not specified
end

if strcmpi(altUnits, 'm') || strcmpi(altUnits, 'meter') || strcmpi(altUnits, 'meters')
    h_m = h;
elseif strcmpi(altUnits, 'ft') || strcmpi(altUnits, 'feet') || strcmpi(altUnits, 'foot')
    h_m = h * ft2m;
else
    error('Altitude units not recognized.')
end

%% Call atmosisa
[T_K, SOS_mps, p_Pa, rho_kg_m3, ~, ~] = atmosisa(h_m);
% Convert to English base units
T_Ra = K2Rankine*T_K; % [Ra]
SOS_fps = mps2fps*SOS_mps; % [ft/s]
p_psf = Pa2psf*p_Pa; % [lbf/ft^2]
rho_lbm_cuft = kg_m3_2lbm_cuft*rho_kg_m3; % [lbm/ft^3]
rho_slug_cuft = rho_lbm_cuft / gc;

%% Density units check
if nargin < 3 || isempty(densityUnits)
    densityUnits = 'slug'; % default to slugs if not specified
end

if strcmpi(densityUnits,'lbm') || strcmpi(densityUnits,'lb') || strcmpi(densityUnits,'pounds') || strcmpi(densityUnits,'lbs')
    rhoEng = rho_lbm_cuft;
elseif strcmpi(densityUnits,'slug') || strcmpi(densityUnits,'slugs') || strcmpi(densityUnits,'sl') || strcmpi(densityUnits,'slug/ft^3')
    rhoEng = rho_slug_cuft;
else
    error('Density units not recognized.')    
end

%% Temperature units check
if nargin < 4 || isempty(tempUnits)
    tempUnits = 'Rankine'; % default to Rankine if not specified
end

if strcmpi(tempUnits,'Ra') || strcmpi(tempUnits,'R') || strcmpi(tempUnits,'Rankine')
    T_Eng = T_Ra;
elseif strcmpi(tempUnits,'F') || strcmpi(tempUnits,'deg F') || strcmpi(tempUnits,'Fahrenheit') || strcmpi(tempUnits,'degrees Fahrenheit')
    T_Eng = T_Ra - 459.67;
else
    error('Temperature units not recognized.')    
end

end

function CD = getCD(CL, M)

CD0_sub = 0.02;
K2_sub = 0.025;

CD0_super = 0.022;
% CD0_super = 0.02;
K2_super = 0.025;

CDR = zeros(size(M));

machCurve = [0        0.85    0.95                               1.05       2             ];
%CD0curve = [0.02     0.02    0.03                               0.045      0.04          ];
CD0curve  = [CD0_sub  CD0_sub 1.125*(CD0_super-CD0_sub)+CD0_sub  CD0_super  0.95*CD0_super];
K2curve   = [K2_sub   K2_sub  1.125*(K2_super-K2_sub)+K2_sub     K2_super   K2_super      ];

CD0 = interp1(machCurve, CD0curve, M, "linear", "extrap");
K2 = interp1(machCurve, K2curve, M, "linear", "extrap");

CD = CD0 + K2 .* (CL.^2) + CDR;

end % end getCD