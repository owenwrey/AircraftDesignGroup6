% V-n diagram 
% based on notes from V-n Diagram Construction Document, 4.3.4.2
clear; clc; close all;

%% Update checklist, MAKE SURE TO DO THIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% things to update include: limit loads, max Mach #, CLmax and min, weight
% fraction, wing loading, thrust-to-weight, MAC, CL_alpha, max sea level
% flight speed and cruise speed. After that, the code should be able to
% auto-update
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Inputs
n_lim_p = 8;            % positive limit load
n_lim_n = -3;           % negative limit load
M = 1.7;                % max Mach number

% aero characteristics 
CL_max = 1.1;           % CL_max
CL_min = -1;            % negative CL_max

beta = .772;            % mid mission weight fraction
WL = 130;               % wing loading, lbf/ft^2
TW = .6;                % thrust to weight
rho_sl = .002377;       % air density, slug/ft^3

% gust line inputs
mac = 13.47;            % mean aerodynamic chord
g = 32.2;               % gravity, ft/s^2
CL_alpha = 4.9;         % from 1/K
U_de = [66,-66,50,-50,25,-25];        % derived gust velocity, ft/s

% max sea-level level flight speed, KEAS
VH = M*sqrt(1.4*1716*518.67)/1.68781;
VD = 1.25*VH;           % design dive speed, KEAS
VC = 573;               % cruise speed 


%% Calculations
% Calc Positive load limits
Vs_pos = sqrt(2*beta*WL/(rho_sl*CL_max));       % positive stall speed
Vs_pos_line = [Vs_pos;Vs_pos];
V = linspace(Vs_pos,VD,100)';                   % create span of V values
q_bar = .5*rho_sl*V.^2;                         % calc q_bar

n_maxaero = zeros(length(V),1);
for i = 1:length(V)
    n_maxaero(i) = q_bar(i)*CL_max/(beta*WL);
    if n_maxaero(i) > n_lim_p
        n_maxaero(i) = n_lim_p;
    end
end

VA = Vs_pos*sqrt(n_lim_p);

% Calc Negative load limits
Vs_neg = sqrt(2*beta*WL/(rho_sl*abs(CL_min)));       % negative stall speed
Vs_neg_line = [Vs_neg;Vs_neg];
V_neg = linspace(Vs_neg,VD,100)';
q_bar = .5*rho_sl*V_neg.^2;

n_minaero = zeros(length(V),1);
for i = 1:length(V)
    n_minaero(i) = q_bar(i)*CL_min/(beta*WL);
    if n_minaero(i) < n_lim_n
        n_minaero(i) = n_lim_n;
    end
end

% Calc gust lines
V_gust = [1;V];

mu_g = 2*WL/(rho_sl*mac*g*CL_alpha);
Kg = .88*mu_g/(5.3 + mu_g);
for i = 1:length(U_de)
    n_gust(:,i) = 1 + (Kg*U_de(i)*V_gust*CL_alpha)/(498*WL);
end


xpts = [Vs_pos;Vs_neg;VA;VH;VD];
labels = {'V_s^+', 'V_s^-', 'VA', 'VH', 'VD'};

%% Plotting
plot(V,n_maxaero,'b', 'LineWidth',2)        % positive load factor
ylim([n_lim_n-1 n_lim_p+1])
ax = gca;
ax.YTick = (n_lim_n-1):1:(n_lim_p+1);
hold on
grid on
xline(0, 'k', 'LineWidth',1)
yline(0, 'k', 'LineWidth',1)
xlabel('Airspeed (KEAS)')
ylabel('Load Factor')
title('Load Factor v. Airspeed')
plot(V_neg,n_minaero,'b', 'LineWidth',2)        % negative load factor
plot(Vs_pos_line,[0;n_maxaero(1)],'b', 'LineWidth',2)   % positive stall line
plot(Vs_neg_line,[0;n_minaero(1)],'b', 'LineWidth',2)   % negtaive stall line
plot([VD;VD],[n_lim_n;n_lim_p],'b', 'LineWidth',2)      % VD line at back
plot([VH;VH],[0;n_lim_p],'k', 'LineStyle','--')         % VA line
plot([VA;VA],[0;n_lim_p],'k', 'LineStyle','--')         % VH line
plot(V_gust,n_gust,'k', 'LineStyle','--')


% create velocity labels 
plot(xpts, zeros(size(xpts)), 'rv', 'MarkerFaceColor','r');
text(xpts(1), 0, labels{1},'VerticalAlignment','bottom', ...
        'HorizontalAlignment','right','FontSize', 10);
for k = 2:numel(xpts)
    text(xpts(k), 0, labels{k}, ...
        'VerticalAlignment','bottom', ...
        'HorizontalAlignment','left', ...
        'FontSize', 10);
end


xpts = [Vs_pos;Vs_neg;VA;VH;VD];
labels = {'V_s^+', 'V_s^-', 'VA', 'VH', 'VD'};

fprintf('The positive load factor stall speed is %1.3f kts \n',Vs_pos)
fprintf('The negative load factor stall speed is %1.3f kts \n',Vs_neg)
fprintf('The maneuver speed is %1.3f kts \n',VA)
fprintf('The maximum design speed is %1.3f kts \n',VH)
fprintf('The maximum dive speed is %1.3f kts \n',VD)