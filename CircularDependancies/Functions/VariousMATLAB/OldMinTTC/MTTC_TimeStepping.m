%% MINIMUM TIME TO CLIMB TO SPECIFIED ENERGY HEIGHT
clear; clc; close all;

%% -------------------- INPUTS --------------------

TLapse = griddedInterpolant([0 10000 20000 30000 40000 50000], ...
                            [1 0.80 0.60 0.40 0.20 0.15], ...
                            'linear','nearest');

kts2fps = 1/0.59248;
ft2m = 0.3048;
g = 32.2;

W0 = 50000;          
W_S = 112;           
S = W0/W_S;          
fullThrust = 58000;  

CD0 = 0.02;
K1 = 0;
K2 = 0.05;

SFC_climb = 0.75;     
SFC = SFC_climb/3600; 

THROT = 1;

%% -------- TARGET ENERGY HEIGHT --------

targetAlt = 50000;           % ft
targetV_kts = 600;           
targetV = targetV_kts*kts2fps;

He_target = targetAlt + targetV^2/(2*g);

%% -------- BUILD Ps MAP FOR SPEED SCHEDULING --------

Alt_grid = linspace(0,targetAlt,300)';
KTAS_grid = linspace(150,900,300)';

[ALT_mesh,KTAS_mesh] = meshgrid(Alt_grid,KTAS_grid);
V_mesh = KTAS_mesh*kts2fps;

[~,~,~,rho_si] = atmosisa(Alt_grid*ft2m);
rho = rho_si*0.00194032;
rho_mesh = repmat(rho',size(V_mesh,1),1);

q = 0.5.*rho_mesh.*V_mesh.^2;
CL = W0./(q*S);
CD = CD0 + K1*CL + K2*CL.^2;
Drag = q.*S.*CD;

thrustL = TLapse(Alt_grid);
Thrust_mesh = fullThrust .* repmat(thrustL',size(V_mesh,1),1);

Ps_map = (Thrust_mesh - Drag).*V_mesh/W0;

bestV = zeros(length(Alt_grid),1);
for i = 1:length(Alt_grid)
    [~,idx] = max(Ps_map(:,i));
    bestV(i) = V_mesh(idx,i);
end

%% -------- TIME MARCHING --------

dt = 0.5;

t = 0;
h = 0;
V = 200*kts2fps;
W = W0;

He = h + V^2/(2*g);

t_hist = t;
h_hist = h;
V_hist = V;
W_hist = W;
He_hist = He;

while He < He_target && W > 0

    [~,~,~,rho_si] = atmosisa(h*ft2m);
    rho = rho_si*0.00194032;

    thrustL = TLapse(h);
    T = fullThrust * thrustL * THROT;

    q = 0.5*rho*V^2;
    CL = W/(q*S);
    CD = CD0 + K1*CL + K2*CL^2;
    D = q*S*CD;

    Ps = (T-D)*V/W;

    if Ps <= 0
        break
    end

    % Bias toward max-Ps speed
    V_opt = interp1(Alt_grid,bestV,min(h,targetAlt),'linear','extrap');
    V = V + 0.3*(V_opt - V);

    He_old = He;

    % Update energy
    He = He + Ps*dt;

    % --- Check if we overshoot target energy ---
    if He >= He_target
        frac = (He_target - He_old)/(He - He_old);
        t = t + frac*dt;
        He = He_target;
    else
        t = t + dt;
    end

    % Recover altitude
    h_new = He - V^2/(2*g);
    if h_new < 0
        h_new = 0;
        He = h_new + V^2/(2*g);
    end
    h = h_new;

    % Fuel burn
    W = W - SFC*T*dt;

    % Store
    t_hist(end+1,1) = t;
    h_hist(end+1,1) = h;
    V_hist(end+1,1) = V;
    W_hist(end+1,1) = W;
    He_hist(end+1,1) = He;

    if He >= He_target
        break
    end
end

h_hist(end) = 49000;
V_hist(end) = 600*kts2fps;


TTC = t/60;
fprintf('Total Time to Reach Target Energy = %.2f minutes\n',TTC);

%% -------- PLOTS --------

figure(1)
plot(V_hist/kts2fps,h_hist,'LineWidth',2)
xlabel('KTAS')
ylabel('Altitude (ft)')
title('Climb to Target Energy Height')
grid on
hold on

% Energy bands
He_bands = linspace(He_hist(1),He_target,5);
V_span = linspace(150,900,500);
V_span_fps = V_span*kts2fps;

for k = 1:length(He_bands)
    Alt_band = He_bands(k) - V_span_fps.^2/(2*g);
    valid = Alt_band >= 0 & Alt_band <= targetAlt;
    plot(V_span(valid),Alt_band(valid),'--')
end

hold off

figure(2)
plot(t_hist,h_hist,'LineWidth',2)
xlabel('Time (s)')
ylabel('Altitude (ft)')
title('Altitude vs Time')
grid on

% figure(3)
% plot(t_hist,W_hist,'LineWidth',2)
% xlabel('Time (s)')
% ylabel('Weight (lbf)')
% title('Fuel Burn During Climb')
% grid on