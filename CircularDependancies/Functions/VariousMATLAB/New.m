% Minimum Time to Climb with fuel burn + Mach
% clear; clc; close all;

% notes: no AB, going 30 kft & 400 kts

%% Inputs

f119perf = f119perfGEN;

% TLapse = griddedInterpolant([0, 10000, 20000, 30000, 40000, 50000], ...
                            % [1, 0.80, 0.60, 0.40, 0.20, 0.15], ...
                            % 'linear','nearest');

rho_SL = 0.0023769;
kts2fps = 1/0.59248;
ft2m = 0.305;

W0 = aircraft.weight.total;             
W = W0;
W_S = aircraft.constants.wingLoading;              
S = aircraft.wing.Area;
% fullThrust = aircraft.engine.thrustMil*2;
targAlt = 40000;
targSpeed = 400;

CD0_sub = 0.02;         
CD0_sup = 0.06;         
K2 = 0.05;
CDR = 0;

% SFC_climb = 0.75/3600;  % lb/lb/s
% THROT = 1;

step = .5;               

%% Storage
dat = zeros(5000,5); % [time, V(kts), h(ft), W(lb), Mach]
dat(1,:) = [0, 150, 0, W0, 0];

i = 1;

while dat(i,3) < targAlt || dat(i,2) < targSpeed
    
    t = dat(i,1);
    V = dat(i,2);
    h = dat(i,3);
    W = dat(i,4);
    
    % Atmosphere
    [~, a, ~, rho] = atmosisa(h*ft2m);
    rho = rho*0.00194032;
    a_fps = a * 3.28084;
    
    % thrust = fullThrust * TLapse(h) * THROT;
    
    % Dynamic speed search
    V_range = linspace(max(100, V-100), V+300, 120);
    
    best_score = -inf;
    V_best = V;
    Ps_best = 0;
    
    % Determine if we should prioritize climb or speed
    if h >= targAlt
        climb_weight = 0; % prioritize speed
    else
        climb_weight = 1; % prioritize climb
    end
    
    for V_try = V_range
        
        V_fps = V_try * kts2fps;
        M = V_fps / a_fps;

        thrust = f119perf.max.fn(h,M)*2;
        
        q = 0.5 * rho * V_fps^2;
        CL = W / (q*S);
        
        % Drag model
        if M < 1
            CD0 = CD0_sub;
        else
            CD0 = CD0_sup;
        end
        
        CD = CD0 + K2*CL^2 + CDR;
        D = q*S*CD;
        
        Ps = (thrust - D)*V_fps / W;
        
        % --- Scoring function ---

        speed_error = V_try - targSpeed;
        speed_penalty = speed_error^2;
        alt_error = targAlt - h;
        blend = max(0, min(1, alt_error / 2000)); % 2000 ft transition band

        % Final score
        score = blend*Ps - (1-blend)*0.01*speed_penalty;
        
        if score > best_score
            best_score = score;
            V_best = V_try;
            Ps_best = Ps;
        end
    end
    
    %% --- State update ---
    
    V_current = V * kts2fps;
    V_target  = V_best * kts2fps;
    
    % Smooth acceleration
    dV = (V_target - V_current)*0.3;
    dV = max(min(dV,50),-50);
    
    % Energy accounting
    dE_dt = Ps_best;
    dE_dV = V_current * dV / 32.2;
    
    dh = max(0, dE_dt - dE_dV);
    
    % Update states
    V_new = V_current + dV*step;
    h_new = max(0, h + dh*step);
    
    % Mach
    M_new = V_new / a_fps;
    
    % Fuel burn
    SFC_climb = f119perf.max.tsfc(h_new,M_new)/3600; 
    thrust = f119perf.max.fn(h,M)*2; 
    W_new = W - SFC_climb * thrust * step;
        
    % Store
    dat(i+1,1) = t + step;
    dat(i+1,2) = V_new / kts2fps;
    dat(i+1,3) = h_new;
    dat(i+1,4) = W_new;
    dat(i+1,5) = M_new;
    
    i = i + 1;
    if i>10e3
        fprintf('shits fucked\n')
        break
    end
end

dat = dat(1:i,:); % trim

%% Plots

figure(2)
plot(dat(:,2), dat(:,3),'b','LineWidth',2)
xlabel('KTAS')
ylabel('Altitude (ft)')
title('Climb Trajectory')
grid on

figure(3)
plot(dat(:,1), dat(:,3),'b','LineWidth',2)
xlabel('Time (s)')
ylabel('Altitude (ft)')
title('Altitude vs Time')
grid on

time = dat(:,1);
speed = dat(:,2);
alt = dat(:,3);

Results = table(time, speed, alt);
% disp(Results)

fprintf('Total time to climb is %.0f minutes \n',(time(end)/60))