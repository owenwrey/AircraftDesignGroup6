clear
close all
clc
%% Script Settings
printStatsAtEachAltitude = 0; % if ==1 your command window gets cluttered

% This what defines what the service ceiling is:
serviceCeilingClimbRate = 100; % [ft/min] (not [fps] here)

% you don't need to change these unless your not reaching the abs ceiling
altRange = 50000; % max altitude tested, if ceiling is below this the script stops
machRangeMax = 2.8; % this is the max speed plotted


%% Aircraft Settings
throtRatio = 1.0; % 1.0 == full afterburner, don't increase past 1

%% Aircraft Parameters
W = 76965; % Aircraft Weight [lb_f]
Sref = 592.04; % Wing Area [ft^2]
Max_Thrust = 44000; % Max Thrust SLS [lb_f]

% Drag Polar
CD0 = 0.02;
K1 = 0;
K2 = 0.129;
CDR = 0;
CLmax = 1.2;

n = 1;

%% Environment Parameters

%[~,~,rho_SL] = atmos_English(0);
[~, ~, ~, rho_SL] = atmosEnglish(0, 'ft', 'slugs', 'R');

altitude = 0:25:altRange; % TO VARY ALTITUDE CHANGE THIS PARAMETER
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
linMach = (0.1:0.0025:machRangeMax);

tbl = table;
tbl.alt = altitude(i)*ones(length(linMach),1);
tbl.Mach = linMach';
tbl.V_TAS = tbl.Mach.*a;
tbl.V_EAS = tbl.V_TAS.*sqrt(sigma);
tbl.KTAS = tbl.V_TAS*fpsToKnots;
tbl.KEAS = tbl.V_EAS*fpsToKnots;
tbl.q = .5.*rho.*tbl.V_TAS.^2;
tbl.CL = W./tbl.q./Sref;
tbl.CD = CD0 + K1.*tbl.CL + K2.*tbl.CL.^2 + CDR;
% tbl.CD = getCD(tbl.CL); % <-- change the name of this function to our drag polar function
tbl.L2D = tbl.CL./tbl.CD;
tbl.CL12_CD = tbl.CL.^(1/2)./tbl.CD;
tbl.Drag = Sref.*tbl.q.*tbl.CD;
tbl.Thrust = (1.11.*(sigma)-0.11).* Max_Thrust .* ones(length(tbl.Mach),1);
% If you want to use the more complex TLapse comment the line above this
% and uncomment the line below
% tbl.Thrust = TLapse(altitude(i), tbl.Mach, throtRatio, 'LBPR').* Max_Thrust .* ones(length(tbl.Mach),1);
tbl.DV = tbl.Drag.*tbl.V_TAS/550;
tbl.TV = tbl.Thrust.*tbl.V_TAS/550;
tbl.hdot = 550*(tbl.TV - tbl.DV)./W; % [fps] (not [fpm] here)
tbl.FPA = real(asind(tbl.hdot./tbl.V_TAS)); % [degree]

% if altitude(i) == 1000 % ft
%     fprintf('i: %u\n', i)
%     fprintf("alt: %f\n", altitude(i))
%     break;
% end


Vx(i) = tbl.V_EAS(max(tbl.FPA) == tbl.FPA)*fpsToKnots; % [KEAS]
Vy(i) = tbl.V_EAS(max(tbl.hdot) == tbl.hdot)*fpsToKnots; % [KEAS]
% find the index where the hdot changes sign:
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
% assign results
VmaxKEAS(i) = tbl.V_EAS(imax) * fpsToKnots;
VmaxKTAS(i) = tbl.V_TAS(imax) * fpsToKnots;
Mmax(i)     = tbl.Mach(imax);

VstallKTAS(i) = sqrt((2/rho) * (W/Sref) *(n/CLmax))*fpsToKnots;
VstallKEAS(i) = VstallKTAS(i) * sqrt(sigma);
Mstall(i)     = VstallKTAS(i)/a;

hdot_x(i) = tbl.hdot(max(tbl.FPA) == tbl.FPA)*60; % {ft/min}
hdot_y(i) = max(tbl.hdot)*60; % {ft/min}
FPA_x(i) = max(tbl.FPA);
FPA_y(i) = tbl.FPA(max(tbl.hdot) == tbl.hdot);

% thrust(i) = tbl.Thrust(1);
% cd(i)     = tbl.CD(1000);

% check if max speed is less than stall speed and assign abs ceil
if VmaxKEAS(i) < VstallKEAS(i)
    fprintf('Absolute Ceiling found (max speed less than stall): %.0f ft\n', altitude(i));
    absCeiling = altitude(i);
    break;
end


%% Print
if printStatsAtEachAltitude == 1
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
fprintf('Service Ceiling: %u ft (%u ft/min climb rate)\n', serviceCeiling, serviceCeilingClimbRate);

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
absCeilingForPlot    = absCeiling * ones(length(altitude),1);
absCeilingVeloForPlot = linspace(0, 1500, length(absCeilingForPlot));

serviceCeilingForPlot = serviceCeiling * ones(length(altitude),1);
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
plot(absCeilingVeloForPlot, absCeilingForPlot, 'LineWidth', 1.75, 'Color', 'k');
plot(serviceCeilingVeloForPlot, serviceCeilingForPlot, 'LineWidth', 1.25, 'LineStyle','--');
grid on;
xlim([100, 1200])
xlabel("Maximum Velocity (KEAS)",Interpreter="latex")
ylabel("Altitude (ft)",Interpreter="latex")
legend('Maximum Speed', 'Stall Speed', 'Absolute Ceiling', 'Service Ceiling');

% MACH vs ALT
figure;
plot(Mmax,altitude,'LineWidth',1.75); hold on;
plot(Mstall, altitude, 'LineWidth',1.75);
plot(absCeilingVeloForPlot, absCeilingForPlot, 'LineWidth', 1.75, 'Color', 'k');
plot(serviceCeilingVeloForPlot, serviceCeilingForPlot, 'LineWidth', 1.25, 'LineStyle','--');
grid on;
xlim([0.15, 1.8]);
xlabel("Maximum Mach Number",Interpreter="latex")
ylabel("Altitude (ft)",Interpreter="latex")
legend('Maximum Speed', 'Stall Speed', 'Absolute Ceiling', 'Service Ceiling');

% KTAS vs ALT
figure;
plot(VmaxKTAS,altitude,'LineWidth',1.75); hold on;
plot(VstallKTAS, altitude, 'LineWidth',1.75);
plot(absCeilingVeloForPlot, absCeilingForPlot, 'LineWidth', 1.75, 'Color', 'k');
plot(serviceCeilingVeloForPlot, serviceCeilingForPlot, 'LineWidth', 1.25, 'LineStyle','--');
grid on;
xlim([100, 1200]);
xlabel("Maximum Velocity (KTAS)",Interpreter="latex")
ylabel("Altitude (ft)",Interpreter="latex")
legend('Maximum Speed', 'Stall Speed', 'Absolute Ceiling', 'Service Ceiling');

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

end