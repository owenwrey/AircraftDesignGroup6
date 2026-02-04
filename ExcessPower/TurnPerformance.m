%% MASTER run file
% NOT FULLY GENERALIZED
clear variables; close all; clc;

%% Unit conversions
fps2kn = 0.592484;
nmi2ft = 6076.12; % naut. miles to feet

% Get atmos conditions

alt = [0, 15e3, 30e3]; % [ft]
% alt = 20e3;
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
alpha_thrust = TLapse(alt); % thrust lapse at 0, 15k, and 30k ft
T = T_SL*alpha_thrust; % row vector of thrusts

% 1 g stall speed
Vs1g = sqrt( (2./(rho)).* (W/Sref)*(1/CLmax) );
% Set Vmin
Vmin = Vs1g;
% Number of velos in array
numV = 50;
%% Find Vmax

%% Find Vmax

for i = 1:num_alt
    
    Vlow = Vs1g(i);
    Vhigh = 6000;

    % if i > 8 
    %     Vhigh = 1700;
    % end
    
    % Excess thrust function (positive when thrust > drag)
    Pexc_func = @(V) T(i) - 0.5 .* rho(i) .* V.^2 .* Sref .* getCD(2.*W./(rho(i).*V.^2.*Sref), V./a(i) );
    
    % Use fzero to find where excess thrust = 0
    options = optimset('TolX', 1e-3);
    try
        Vmax(i) = fzero(Pexc_func, [Vlow, Vhigh], options);
    catch ME
        warning('fzero failed for altitude index %d: %s', i, ME.message);
        % Fallback: use upper bound or NaN
        Vmax(i) = NaN;
    end
    
end


for i = 1:num_alt
    V(i,:) = linspace(Vs1g(i), Vmax(i), numV); % ignore error
end
% get q matrix
for i = 1:num_alt
    q(i,:) = 0.5*rho(i)*V(i,:).^2;
end

%% Compute Loads

% max load factor before stall
for i = 1:num_alt
    for j = 1:numV
        nmaxInst(i,j) = q(i,j).*CLmax./(W/Sref);

        if nmaxInst(i,j) > nmaxStruct
            nmaxInst(i,j) = nmaxStruct;
        end
    end
end


for i = 1:num_alt
    for j = 1:numV
        qi = q(i,j);
        nLow = 1;
        nHigh = nmaxInst(i,j);
        if nHigh <= nLow
            nmaxSus(i,j) = NaN;
            continue;
        end
        
        for k = 1:50
            nMid = (nHigh+nLow)/2;
            CLi = nMid*W/(qi*Sref);

            if CLi > CLmax
                nHigh = nMid;
                continue;
            end
            CD = getCD(CLi, V(j)/a(i));
            D = qi*Sref*CD;
            % D = qi*Sref*(CD0 + K1*CLi + K2*CLi^2);
            
            if T(i) - D > 0
                nLow = nMid;
            else
                nHigh = nMid;
            end
        end % end for k
        nmaxSus(i,j) = nMid;
        nmaxSus(:,1) = 1.00; % nmax at stall must be 1
    end % end for j
end % end for i

%% Compute other factors

for i = 1:num_alt
    nmaxSus(i,numV) = floor(nmaxSus(i,numV));
    for j = 1:numV
        phiSus(i,j) = acos(1/nmaxSus(i,j))*180/pi;
        phiInst(i,j) = real(acos(1/nmaxInst(i,j))*180/pi); % small imaginary component for some reason

        R(i,j) = V(i,j)^2/(gc*sqrt(nmaxSus(i,j)^2-1))/nmi2ft;

        psiDotSus(i,j) =  gc*sqrt(nmaxSus(i,j)^2-1)/V(i,j)*180/pi;
        psiDotInst(i,j) = real(gc*sqrt(nmaxInst(i,j)^2-1)/V(i,j))*180/pi;
    end
end

%% Plot
% convert to knots
V = fps2kn*V;

% figure('Units','normalized','Position',[0.05 0.05 0.9 0.85]);
figure;
%{
colors = 0.8 * [
    0.0 0.0 1.0   % blue
    0.0 0.3 1.0
    0.0 0.6 1.0
    0.0 0.9 1.0
    0.0 1.0 0.5   % cyan-green
    0.0 1.0 0.0   % green
    0.5 1.0 0.0
    1.0 1.0 0.0   % yellow
    1.0 0.5 0.0
    1.0 0.0 0.0   % red
];
%}
colors = 0.85* [
    0.0 0.0 1.0
    0.0 1.0 0.0
    1.0 0.0 0.0];

lw = 1.6;

for i = 1:num_alt
    % load factor (n)
    subplot(4,1,1);
    plot(V(i,:),nmaxSus(i,:) ,'Color',colors(i,:),'LineWidth',lw); hold on;
    plot(V(i,:),nmaxInst(i,:),'Color',colors(i,:), 'LineStyle','--','LineWidth',lw);
    ylim([0 9]);
    % xlim([0 1400]);
    xlabel('KTAS');
    ylabel('Load Factor');

    if i == num_alt
        % legend('0 ft Sustained','0 ft Instantaneous', ...
        %        '5k ft Sustained','5k ft Instantaneous', ...
        %        '10k ft Sustained','10k ft Instantaneous', ...
        %        '15k ft Sustained','15k ft Instantaneous', ...
        %        '20k ft Sustained','20k ft Instantaneous', ...
        %        '25k ft Sustained','25k ft Instantaneous', ...
        %        '30k ft Sustained','30k ft Instantaneous', ...
        %        '35k ft Sustained','35k ft Instantaneous', ...
        %        '40k ft Sustained','40k ft Instantaneous', ...
        %        '45k ft Sustained','45k ft Instantaneous');
        legend('0 ft Sustained',    '0 ft Instantaneous', ...
               '15k ft Sustained','15k ft Instantaneous', ...
               '30k ft Sustained','30k ft Instantaneous');
    end

    % % phi (bank angle)
    % subplot(4,1,2);
    % plot(V(i,:),phiSus(i,:) ,'Color',colors(i,:)); hold on;
    % plot(V(i,:),phiInst(i,:),'Color',colors(i,:), 'LineStyle','--');
    % xlabel('KTAS');
    % ylabel('Bank Angle (deg)');

    % psiDot (turn rate)
    subplot(4,1,2);
    plot(V(i,:),psiDotSus(i,:) ,'Color',colors(i,:),'LineWidth',lw); hold on;
    plot(V(i,:),psiDotInst(i,:),'Color',colors(i,:), 'LineStyle','--','LineWidth',lw);
    xlabel('KTAS');
    ylabel('Turn Rate (deg/s)');

    % turn radius (nmi)
    subplot(4,1,3);
    plot(V(i,:),R(i,:) ,'Color',colors(i,:),'LineWidth',lw); hold on;
    xlabel('KTAS');
    ylabel('Turn Radius (NMI)');
    % ylim([0,30]);

end


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