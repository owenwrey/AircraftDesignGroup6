%% MASTER run file
% NOT FULLY GENERALIZED
clear variables; close all; clc;

%% Unit conversions
fps2kn = 0.592484;
nmi2ft = 6076.12; % naut. miles to feet

% Get atmos conditions

alt = 0:5000:50000; % [ft] - general altitude array
[~, a, ~, rho] = atmosEnglish(alt, 'ft', 'slug');
num_alt = length(alt);

gc = 32.174; % [lbm/slug]

%% Aircraft Characteristics
% W0 = 75250; % [lbf]
% beta = 0.8485;
% W = W0*beta;
% Sref = 578.846; % [ft^2]

W0 = 84127; % Aircraft Weight [lb_f]
Sref = W0/115; % Wing Area [ft^2]
W = 72614;
% W = W0;
T_SL = 58000; % Max Thrust SLS [lb_f]

CLmax = 1.2;
% CD0 = 0.02;
% K1 = 0;
% K2 = 0.025;

nmaxStruct = 8.0; % [g]

%% Thrust array

TLapse = griddedInterpolant([0, 10000, 20000, 30000, 40000, 50000], [1, 0.80, 0.60, ...
0.40, 0.20, 0.15],'linear','nearest');
alpha_thrust = TLapse(alt); % thrust lapse at altitudes
T = T_SL*alpha_thrust; % row vector of thrusts

% 1 g stall speed
Vs1g = sqrt( (2./(rho)).* (W/Sref)*(1/CLmax) );
% Set Vmin
Vmin = Vs1g;
% Number of velos in array
numV = 50;

%% Find Vmax and compute velocity/dynamic pressure arrays
Vmax = zeros(1, num_alt);
V = zeros(num_alt, numV);
q = zeros(num_alt, numV);

for i = 1:num_alt
    Vlow = Vs1g(i);
    Vhigh = 6000;
    
    % Excess thrust function (positive when thrust > drag)
    Pexc_func = @(V) T(i) - 0.5 .* rho(i) .* V.^2 .* Sref .* getCD(2.*W./(rho(i).*V.^2.*Sref), V./a(i) );
    
    % Use fzero to find where excess thrust = 0
    options = optimset('TolX', 1e-3);
    try
        Vmax(i) = fzero(Pexc_func, [Vlow, Vhigh], options);
    catch ME
        warning('fzero failed for altitude index %d: %s', i, ME.message);
        Vmax(i) = NaN;
    end
    
    % Create velocity array and compute dynamic pressure
    V(i,:) = linspace(Vs1g(i), Vmax(i), numV);
    q(i,:) = 0.5*rho(i)*V(i,:).^2;
end

%% Compute Loads
% Preallocate arrays
nmaxInst = zeros(num_alt, numV);
nmaxSus = zeros(num_alt, numV);
phiSus = zeros(num_alt, numV);
phiInst = zeros(num_alt, numV);
R = zeros(num_alt, numV);
psiDotSus = zeros(num_alt, numV);
psiDotInst = zeros(num_alt, numV);

% Compute instantaneous and sustained load factors
for i = 1:num_alt
    % Vectorized instantaneous load factor calculation
    nmaxInst(i,:) = min(q(i,:) .* CLmax ./ (W/Sref), nmaxStruct);
    
    % Sustained load factor calculation
    for j = 1:numV
        qi = q(i,j);
        nLow = 1;
        nHigh = nmaxInst(i,j);
        
        if nHigh <= nLow
            nmaxSus(i,j) = NaN;
            continue;
        end
        
        % Function to find sustained load factor
        % At sustained n, thrust must equal drag: T = D
        % D = q*Sref*CD, where CL = n*W/(q*Sref)
        thrust_minus_drag = @(n) T(i) - qi*Sref*getCD(n*W/(qi*Sref), V(i,j)/a(i));
        
        % Check if there's a sign change (zero crossing exists)
        f_low = thrust_minus_drag(nLow);
        f_high = thrust_minus_drag(nHigh);
        
        if f_low * f_high > 0
            % No sign change - aircraft can sustain maximum instantaneous n
            % (thrust > drag across entire range)
            nmaxSus(i,j) = nHigh;
        else
            % Sign change exists - use fzero to find the crossing
            options = optimset('TolX', 1e-4);
            nmaxSus(i,j) = fzero(thrust_minus_drag, [nLow, nHigh], options);
        end
    end
end

% Enforce nmax = 1 at stall speed (first velocity point)
nmaxSus(:,1) = 1.00;

%% Compute other factors (bank angle, turn rate, turn radius)
for i = 1:num_alt
    nmaxSus(i,numV) = floor(nmaxSus(i,numV));
    
    for j = 1:numV
        % Bank angles
        phiSus(i,j) = acos(1/nmaxSus(i,j))*180/pi;
        phiInst(i,j) = real(acos(1/nmaxInst(i,j))*180/pi);
        
        % Turn radius (nmi)
        R(i,j) = V(i,j)^2/(gc*sqrt(nmaxSus(i,j)^2-1))/nmi2ft;
        
        % Turn rates (deg/s)
        psiDotSus(i,j) = gc*sqrt(nmaxSus(i,j)^2-1)/V(i,j)*180/pi;
        psiDotInst(i,j) = real(gc*sqrt(nmaxInst(i,j)^2-1)/V(i,j))*180/pi;
    end
end

%% Plot
% convert to knots
V_kts = fps2kn*V;

figure;
hold on;
grid on;
box on;

% Color map for altitudes (blue to red gradient)
colors = jet(num_alt);

% Plot sustained turn rate curves for each altitude
for i = 1:num_alt
    plot(V_kts(i,:), psiDotSus(i,:), 'Color', colors(i,:), 'LineWidth', 2, ...
        'DisplayName', sprintf('%d ft', alt(i)));
end

% Add iso-load factor lines (constant n lines)
n_iso = [2, 4, 6, 8]; % Load factors to show
V_iso = linspace(100, 1200, 100); % Velocity range for iso lines

for n = n_iso
    psiDot_iso = gc*sqrt(n^2-1)./(V_iso/fps2kn)*180/pi;
    plot(V_iso, psiDot_iso, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    % Add text label
    text(V_iso(end-10), psiDot_iso(end-10), sprintf('n = %d', n), ...
        'FontSize', 9, 'BackgroundColor', 'white', 'EdgeColor', 'k');
end

% Add iso-turn radius lines (constant R lines)
R_iso = [0.25, 0.5, 1, 2, 5, 10, 20]; % Turn radii in nmi to show

for R_nmi = R_iso
    R_ft = R_nmi * nmi2ft;
    psiDot_iso = V_iso/fps2kn / R_ft * 180/pi;
    plot(V_iso, psiDot_iso, 'k:', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    % Add text label
    text(V_iso(end-5), psiDot_iso(end-5), sprintf('R = %d nmi', R_nmi), ...
        'FontSize', 9, 'BackgroundColor', 'white', 'EdgeColor', 'k');
end

xlabel('Velocity (KTAS)', 'FontSize', 12);
ylabel('Turn Rate (deg/s)', 'FontSize', 12);
title('Sustained Turn Performance', 'FontSize', 14);
legend('Location', 'northeast');
xlim([0 1200]);
ylim([0 25]);

hold off;


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