%% Script Settings
printStatsAtEachAltitude = false; % if ==1 your command window gets cluttered

% This what defines what the service ceiling is:
serviceCeilingClimbRate = 500; % [ft/min] (not [fps] here)

% you don't need to change these unless you're not reaching the abs ceiling
altRange = 65000; % max altitude tested, if ceiling is below this the script stops
altStep = 250; % [ft]
machRangeMax = 3.0; % this is the max speed plotted
machStep = 0.0025;


%% Aircraft Settings
throtRatio = 1.0; % 1.0 == full afterburner, don't increase past 1

%% Aircraft Parameters
% W0 = aircraft.weight.midMission; % Aircraft Weight [lb_f]
W0 = 56600; % Aircraft Weight [lb_f]

% W_S = aircraft.constants.wingLoading;
W_S = 102;

Sref = W0/W_S; % Wing Area [ft^2]

% beta = aircraft.TimeStepTable.beta(floor(end/2))
% beta = 0.86;
beta = 1;

W = beta*W0;
% W = W0;

% Max_Thrust = 2*aircraft.engine.thrust;
Max_Thrust = 73000; % Max Thrust SLS [lb_f]

% f119perf.max.fn = griddedInterpolant(ALTq',Mnq',Thrustq');
% f119perf.max.tsfc = griddedInterpolant(ALTq',Mnq',TSFCq');
engPerf = f119perfGEN;

% Drag Polar
% CD0 = 0.02;
% K1 = 0;
% K2 = 0.129;
CDR = 0;
CLmax = 1.2;

n = 1;

%% Environment Parameters

%[~,~,rho_SL] = atmos_English(0);
[~, ~, ~, rho_SL] = atmosEnglish(0, 'ft', 'slugs', 'R');

altitude = 0:altStep:altRange; % TO VARY ALTITUDE CHANGE THIS PARAMETER
% you cant set this over your service ceiling

fpsToKnots = 0.592484;

%% Loop Initialization

Vx       = zeros(1,length(altitude));
Vy       = zeros(1,length(altitude));

VmaxKEAS = zeros(1,length(altitude));
VmaxKTAS = zeros(1,length(altitude));
Mmax     = zeros(1,length(altitude));

VstallKTAS = zeros(1,length(altitude));
VstallKEAS = zeros(1,length(altitude));
Mstall     = zeros(1,length(altitude));

hdot_x   = zeros(1,length(altitude));
hdot_y   = zeros(1,length(altitude));
FPA_x    = zeros(1,length(altitude));
FPA_y    = zeros(1,length(altitude));

Tcell = cell(length(altitude),1);


for i = 1:length(altitude)
%[T,~,rho] = atmos_English(altitude(i));
[T, a, ~, rho] = atmosEnglish(altitude(i), 'ft', 'slugs', 'R');
sigma = rho/rho_SL;
%a = sqrt(1.4.*1716.*T);

%% Acceleration Analysis
linMach = (0.2:machStep:machRangeMax);

tbl = table;
tbl.alt = altitude(i)*ones(length(linMach),1);
tbl.Mach = linMach';
tbl.V_TAS = tbl.Mach.*a;
tbl.V_EAS = tbl.V_TAS.*sqrt(sigma);
tbl.KTAS = tbl.V_TAS*fpsToKnots;
tbl.KEAS = tbl.V_EAS*fpsToKnots;
tbl.q = .5.*rho.*tbl.V_TAS.^2;
tbl.CL = W./tbl.q./Sref;
tbl.CD0 = aircraft.aero.cd0_strike_interp(altitude(i) .* ones(size(tbl.Mach)), tbl.Mach); % assign CD0s from all the Machs
% ff
for j = 2:length(tbl.CD0)
    if tbl.CD0(j) > 0.06
        % tbl.CD0(j) = tbl.CD0(j-1) - 0.00001;
        tbl.CD0(j) = 0.060;
    end
end

tbl.K2 = getK2(tbl.Mach);
tbl.K2 = 0.0*ones(length(linMach),1);
tbl.CD = tbl.CD0 + tbl.K2.*tbl.CL.^2 + CDR;
% tbl.CD = getCD(tbl.CL, tbl.Mach); % <-- change the name of this function to our drag polar function
tbl.L2D = tbl.CL./tbl.CD;
tbl.CL12_CD = tbl.CL.^(1/2)./tbl.CD;
tbl.Drag = Sref.*tbl.q.*tbl.CD;
tbl.Thrust = (1.11.*(sigma)-0.11).* Max_Thrust .* ones(length(tbl.Mach),1);
% tbl.Thrust = engPerf.max.fn(altitude(i).* ones(size(tbl.Mach)), tbl.Mach);
% If you want to use the more complex TLapse comment the line above this
% and uncomment the line below
% tbl.Thrust = TLapse(altitude(i), tbl.Mach, throtRatio, 'LBPR').* Max_Thrust .* ones(length(tbl.Mach),1);
tbl.DV = tbl.Drag.*tbl.V_TAS/550; % [hp]
tbl.TV = tbl.Thrust.*tbl.V_TAS/550; % [hp]
tbl.hdot = 550*(tbl.TV - tbl.DV)./W; % [fps] (not [fpm] here)
% tbl.FPA = real(asind(tbl.hdot./tbl.V_TAS)); % [degree]
tbl.FPA = asind(tbl.hdot./tbl.V_TAS); % [degree]

% if altitude(i) == 1000 % ft
%     fprintf('i: %u\n', i)
%     fprintf("alt: %f\n", altitude(i))
%     break;
% end


% take highest V_EAS if there is a tie for 
Vx(i) = tbl.V_EAS(find(max(tbl.FPA) == tbl.FPA, 1, 'last'))*fpsToKnots; % [KEAS]
Vy(i) = tbl.V_EAS(find(max(tbl.hdot) == tbl.hdot, 1, 'last'))*fpsToKnots; % [KEAS]

% find the index where hdot changes sign:
idx = find(tbl.hdot(2:end) .* tbl.hdot(1:end-1) < 0) + 1; 
% find the row where V_EAS is largest among sign-change points (there's probably only two options)
[~, k] = max(tbl.V_EAS(idx));
% index of speed where you can't climb anymore (hdot ≈0 KEAS) == max speed
imax = idx(k);

if isempty(imax)
    fprintf('Absolute Ceiling: %.0f ft\n', altitude(i-1));
    absCeiling = altitude(i-1);
    break;
end

% if ~isempty(imax)
%% Assign results
VmaxKEAS(i) = tbl.V_EAS(imax) * fpsToKnots;
VmaxKTAS(i) = tbl.V_TAS(imax) * fpsToKnots;
Mmax(i)     = tbl.Mach(imax);
% end
VstallKTAS(i) = sqrt((2/rho) * (W/Sref) *(n/CLmax))*fpsToKnots;
VstallKEAS(i) = VstallKTAS(i) * sqrt(sigma);
Mstall(i)     = (VstallKTAS(i)/fpsToKnots)/a;

hdot_x(i) = tbl.hdot(max(tbl.FPA) == tbl.FPA)*60; % {ft/min}
hdot_y(i) = max(tbl.hdot)*60; % {ft/min}
FPA_x(i) = max(tbl.FPA);
FPA_y(i) = tbl.FPA(max(tbl.hdot) == tbl.hdot);

% thrust(i) = tbl.Thrust(1);
% cd(i)     = tbl.CD(555);

% check if max speed is less than stall speed and assign abs ceil
if VmaxKEAS(i) < VstallKEAS(i)
    fprintf('Absolute Ceiling found (max speed less than stall): %.0f ft\n', altitude(i));
    absCeiling = altitude(i);
    break;
end


%% Print
if printStatsAtEachAltitude % && (altitude(i) == 0 || altitude(i) == 10000 || altitude(i) == 20000 || altitude(i) == 30000 || altitude(i) == 40000 || altitude(i) == 50000)
fprintf('============================\n')
fprintf('At altitude %.0f ft:\n',altitude(i))

fprintf('===== Best Rate Climb =====\n')
fprintf('Speed (Vy): %f KEAS\n',Vy(i))
fprintf('Rate: %f ft/min\n',hdot_y(i))
fprintf('Angle: %f deg\n',FPA_y(i))

fprintf('===== Best Angle Climb =====\n')
fprintf('Speed (Vx): %f KEAS\n',Vx(i))
fprintf('Rate: %f ft/min\n',hdot_x(i))
fprintf('Angle: %f deg\n\n',FPA_x(i))
end % end if

Tcell{i} = tbl;

end % "for i = 1:length(altitude)" end
serviceCeiling = altitude(find(hdot_y >= serviceCeilingClimbRate, 1, 'last')); % find highest altitud where hdot_y is still 200 ft/min
fprintf('Service Ceiling:  %u ft (%u ft/min climb rate)\n', serviceCeiling, serviceCeilingClimbRate);
fprintf('Maximum velocity: %.0f KEAS\n', max(VmaxKEAS))
fprintf('Maximum velocity: %.0f KTAS\n', max(VmaxKTAS))
fprintf('Maximum Mach Num: %.2f\n', max(Mmax))

%% Prep data for plots

altitude = altitude(1:i);

Vx       = Vx(1:i);
Vy       = Vy(1:i);

VmaxKEAS = VmaxKEAS(1:i);
VmaxKTAS = VmaxKTAS(1:i);
Mmax     = Mmax(1:i);

VstallKTAS = VstallKTAS(1:i);
VstallKEAS = VstallKEAS(1:i);
Mstall     = Mstall(1:i);

hdot_x   = hdot_x(1:i);
hdot_y   = hdot_y(1:i);
FPA_x    = FPA_x(1:i);
FPA_y    = FPA_y(1:i);


% check if absolute ceiling was reached
if ~exist("absCeiling", "var")
    absCeiling = altRange;
    warning('Max altitude checked is not the absolute ceiling');
end

% Trim to valid altitudes (below absolute ceiling)
validIdx = altitude < absCeiling;

altitude   = altitude(validIdx);
Vx         = Vx(validIdx);
Vy         = Vy(validIdx);
VmaxKEAS   = VmaxKEAS(validIdx);
VmaxKTAS   = VmaxKTAS(validIdx);
Mmax       = Mmax(validIdx);
VstallKTAS = VstallKTAS(validIdx);
VstallKEAS = VstallKEAS(validIdx);
Mstall     = Mstall(validIdx);
hdot_x     = hdot_x(validIdx);
hdot_y     = hdot_y(validIdx);
FPA_x      = FPA_x(validIdx);
FPA_y      = FPA_y(validIdx);

absCeilingForPlot    = absCeiling * ones(length(altitude),1)';
absCeilingVeloForPlot = linspace(0, 1500, length(absCeilingForPlot));

serviceCeilingForPlot = serviceCeiling * ones(length(altitude),1)';
serviceCeilingVeloForPlot = linspace(0, 1500, length(serviceCeilingForPlot));


%% Plot
% % KTAS vs. CLIMB RATE
% figure;
% plot(tbl.V_TAS*fpsToKnots,tbl.hdot*60, 'LineWidth',1.75)
% grid on
% xlabel("Maximum Velocity (KTAS)",Interpreter="latex")
% ylabel("Climb Rate ($\frac{ft}{min}$)",Interpreter="latex")

% KEAS vs ALT
figure;
plot(VmaxKEAS,altitude,'LineWidth',1.75); hold on;
plot(VstallKEAS,altitude,'LineWidth',1.75)
yline(absCeiling, 'k-', 'LineWidth', 1.75)
yline(serviceCeiling, '--', 'LineWidth', 1.25, 'Color', [0 0.5 0]);
grid on;
xlim([100, max(VmaxKEAS)*1.075])
ylim([0 absCeiling*1.1])
ax = gca;
ax.YAxis.Exponent = 3; % 1000s feet
xlabel("Maximum Velocity (KEAS)",Interpreter="latex")
ylabel("Altitude (ft)",Interpreter="latex")
legend('Maximum Speed', 'Stall Speed', 'Absolute Ceiling', 'Service Ceiling');

% MACH vs ALT
figure;
plot(Mmax,altitude,'LineWidth',1.75); hold on;
plot(Mstall, altitude, 'LineWidth',1.75);
yline(absCeiling, 'k-', 'LineWidth', 1.75)
yline(serviceCeiling, '--', 'LineWidth', 1.25, 'Color', [0 0.5 0]);
grid on
xlim([0.15, max(Mmax)*1.075])
ylim([0 absCeiling*1.1])
ax = gca;
ax.YAxis.Exponent = 3; % 1000s feet
xlabel("Maximum Mach Number",Interpreter="latex")
ylabel("Altitude (ft)",Interpreter="latex")
legend('Maximum Speed', 'Stall Speed', 'Absolute Ceiling', 'Service Ceiling');

% KTAS vs ALT
figure;
plot(VmaxKTAS,altitude,'LineWidth',1.75); hold on;
plot(VstallKTAS, altitude, 'LineWidth',1.75);
yline(absCeiling, 'k-', 'LineWidth', 1.75)
yline(serviceCeiling, '--', 'LineWidth', 1.25, 'Color', [0 0.5 0]);
grid on;
xlim([100, max(VmaxKTAS)*1.075])
ylim([0 absCeiling*1.1])
ax = gca;
ax.YAxis.Exponent = 3; % 1000s feet
xlabel("Maximum Velocity (KTAS)",Interpreter="latex")
ylabel("Altitude (ft)",Interpreter="latex")
legend('Maximum Speed', 'Stall Speed', 'Absolute Ceiling', 'Service Ceiling');

% % Hdot_y vs ALT
% figure;
% plot(hdot_y,altitude,'LineWidth',1.75); hold on;
% grid on;
% xlabel("Climb Rate ($\frac{ft}{min}$)",Interpreter="latex")
% ylabel("Altitude (ft)",Interpreter="latex")

% figure;
% plot(thrust, altitude(1:end-1))
% xlabel('thrust')
% ylabel('altitude')
% 
% figure;
% plot(cd, altitude(1:end-1))
% xlabel('cd')
% ylabel('altitude')


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

% Call atmosisa
[T_K, SOS_mps, p_Pa, rho_kg_m3, ~, ~] = atmosisa(h_m);

% Convert to English base units
T_Ra = K2Rankine*T_K; % [Ra]
SOS_fps = mps2fps*SOS_mps; % [ft/s]
p_psf = Pa2psf*p_Pa; % [lbf/ft^2]
rho_lbm_cuft = kg_m3_2lbm_cuft*rho_kg_m3; % [lbm/ft^3]
rho_slug_cuft = rho_lbm_cuft / gc;

% Density units check
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

% Temperature units check
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

end % end atmosEnglish

function K2 = getK2(M)

% CD0_sub = 0.02;
K2_sub = 0.025;

% CD0_super = 0.022;
K2_super = 0.025;

% CDR = zeros(size(M));

machCurve = [0        0.85    0.95                               1.05       2             ];
%CD0curve = [0.02     0.02    0.03                               0.045      0.04          ];
% CD0curve  = [CD0_sub  CD0_sub 1.125*(CD0_super-CD0_sub)+CD0_sub  CD0_super  0.95*CD0_super];
K2curve   = [K2_sub   K2_sub  1.125*(K2_super-K2_sub)+K2_sub     K2_super   K2_super      ];

% CD0 = interp1(machCurve, CD0curve, M, "linear", "extrap");
K2 = interp1(machCurve, K2curve, M, "linear", "extrap");

% CD = CD0 + K2 .* (CL.^2) + CDR;

end % end getCD

