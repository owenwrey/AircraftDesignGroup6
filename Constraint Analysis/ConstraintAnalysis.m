clear % testing testing
clc
close all
%% Aircraft Setup

%%%%%%%% Aerodynamic Parameters %%%%%%%%

CD_0 = 0.02; % zero-lift drag coeficient
k1 = 0; % drag polar constant 1
k2 = .2; % drag polar constant 2
CD_R = 0; % resultant drag coeficient

TOP = 300; % Take Off Parameter

%% Master Equation Constraint Analysis

W2S = (10:10:350).*47.88025888888984; % Wing Loadings {N/m^2}
[~,~,~,rho_SL] = atmosisa(0,"extended","on"); % Air Density {kg/m^3}

%%%%%%%% Cruise %%%%%%%%

alt = 0 * 0.3048; % Altitude {m}
n = 1; % Load Factor
M = 1.7; % Mach number
[~,a,~,rho] = atmosisa(alt,"extended","on"); % Air Density {kg/m^3}
v = M.*a; % Velocity {m/s}
ddt_h = 0; % Climb Speed {m/s}
ddt_v = 0; % Acceleration {m/s^2}
alpha = 1.11*(rho/rho_SL)-0.11; % Thrust Lapse
beta = 1; % Weight Lapse

T2W_Cruise = Master_Eqn(CD_0,k1,k2,CD_R,n,v,rho,ddt_h,ddt_v,alpha,beta,W2S);

%%%%%%%% Climb %%%%%%%%

alt = 0; % Altitude {m}
n = 1; % Load Factor
M = 0.8; % Mach number
[~,a,~,rho] = atmosisa(alt,"extended","on"); % Air Density {kg/m^3}
v = M.*a; % Velocity {m/s}
ddt_h = 28000*0.3048/60; % Climb Speed {m/s}
ddt_v = 0; % Acceleration {m/s^2}
alpha = 1.11*(rho/rho_SL)-0.11; % Thrust Lapse
beta = 1.0; % Weight Lapse

T2W_Climb = Master_Eqn(CD_0,k1,k2,CD_R,n,v,rho,ddt_h,ddt_v,alpha,beta,W2S);

%%%%%%% Ceiling %%%%%%%%

alt = 46250*0.3048; % Altitude {m}
n = 1; % Load Factor
M = 1.4; % Mach number
[T,a,~,rho] = atmosisa(alt,"extended","on"); % Air Density {kg/m^3}
v = M.*a; % Velocity {m/s}
ddt_h = 500*0.3048/60; % Climb Speed {m/s}
ddt_v = 0; % Acceleration {m/s^2}
alpha = 1.11*(rho/rho_SL)-0.11; % Thrust Lapse
beta = 1; % Weight Lapse

T2W_Ceiling = Master_Eqn(CD_0,k1,k2,CD_R,n,v,rho,ddt_h,ddt_v,alpha,beta,W2S);

%%%%%%% Sustained Turn %%%%%%%%

alt = 20000*0.3048; % Altitude {m}
n = 5.346; % Load Factor % 10 deg sustained turn
M = 0.9; % Mach number
[~,a,~,rho] = atmosisa(alt,"extended","on"); % Air Density {kg/m^3}
v = M.*a; % Velocity {m/s}
ddt_h = 0; % Climb Speed {m/s}
ddt_v = 0; % Acceleration {m/s^2}
alpha = 1.11*(rho/rho_SL)-0.11; % Thrust Lapse
beta = 1; % Weight Lapse

T2W_SustainedTurn = Master_Eqn(CD_0,k1,k2,CD_R,n,v,rho,ddt_h,ddt_v,alpha,beta,W2S);

%%%%%%% Max Speed %%%%%%%%

alt = 30000*0.3048; % Altitude {m}
n = 1; % Load Factor
M = 2; % Mach number
[~,a,~,rho] = atmosisa(alt,"extended","on"); % Air Density {kg/m^3}
v = M.*a; % Velocity {m/s}
ddt_h = 0; % Climb Speed {m/s}
ddt_v = 0; % Acceleration {m/s^2}
alpha = 1.11*(rho/rho_SL)-0.11; % Thrust Lapse
beta = 1; % Weight Lapse

T2W_MaxSpeed = Master_Eqn(CD_0,k1,k2,CD_R,n,v,rho,ddt_h,ddt_v,alpha,beta,W2S);

%% Takeoff Constraint Analysis %%

%%%%%%% TakeOff %%%%%%%%

alt = 0; % Altitude {m}
[~,~,rho] = atmosisa(alt); % Air Density {kg/m^3}
Cl = 1.6/cosd(24);

T2W_Takeoff = TakeoffConstraintAnalysis(TOP,rho/1.225,Cl,W2S);

%% Vertical Constraint Analysis %%

%%%%%%% Landing %%%%%%%%

alt = 0; % Altitude {m}
n = 1; % Load Factor
beta = 1; % Weight Lapse
Cl = 2.0/cosd(24); % Lift Coefficient
k = 1.2; % Safety Factor
v = 188*1.1*0.514444; % Velocity {m/s}
[~,~,~,rho] = atmosisa(alt,"extended","on"); % Air Density {kg/m^3}

W2S_Landing = verticalConstraintAnalysis(n,beta,Cl,k,v,rho);

%%%%%%% Stall %%%%%%%%

alt = 0; % Altitude {m}
n = 1; % Load Factor
beta = 1; % Weight Lapse
Cl = 1.1; % Lift Coefficient
k = 1.0; % Safety Factor
v = 188*0.514444; % Velocity {m/s}
[~,~,~,rho] = atmosisa(alt,"extended","on"); % Air Density {kg/m^3}

W2S_Stall = verticalConstraintAnalysis(n,beta,Cl,k,v,rho);

%% Plotting
close all
hold on
  ax = gca;
ax.ColorOrder = lines(13);
set(ax,'DefaultLineLineWidth', 1.2);
box on
hold on
  
plot(W2S.*0.02088547,T2W_Cruise,'LineStyle','-')
plot(W2S.*0.02088547,T2W_Climb,'LineStyle','-')
plot(W2S.*0.02088547,T2W_Ceiling,'LineStyle','-')
plot(W2S.*0.02088547,T2W_SustainedTurn,'LineStyle','-')
plot(W2S.*0.02088547,T2W_MaxSpeed,'LineStyle','-')
plot(W2S.*0.02088547,T2W_Takeoff,'LineStyle','-')
xl1=xline(W2S_Landing.*0.02088547,'-r');
xl2=xline(W2S_Stall.*0.02088547,'-b');
plot(125, 1.3,"o")
xlabel("Wing Loading (lb/ft^2)")
ylabel("Thrust-to-Weight")
ylim([0,2.5])

legend(["Cruise","Climb","Celing","Turn","Max Speed","Takeoff", "Landing", "Stall", "Point"])

%% COMMENT %%
