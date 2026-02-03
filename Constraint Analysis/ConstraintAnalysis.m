clear
clc
close all

%% Aircraft Setup

%%%%%%%% Aerodynamic Parameters %%%%%%%%

CD_0 = 0.02; 
k1 = 0; 
k2 = .05; 
CD_R = 0; 

TOP = 300; 

%% ISA + Delta T Setup

DeltaT = 20;           % ISA deviation [K]
gamma = 1.4;
R = 287.05;

%% Master Equation Constraint Analysis

W2S = (10:10:350).*47.88025888888984; % {N/m^2}
[~,~,~,rho_SL] = atmosisa(0,"extended","on"); % ISA sea-level density

%%%%%%%% Cruise %%%%%%%%

alt = 30000 * 0.3048;
n = 1;
M = 1.2;

[T_ISA,a_ISA,~,rho_ISA] = atmosisa(alt,"extended","on");
T = T_ISA + DeltaT;
a = sqrt(gamma*R*T);
rho = rho_ISA .* (T_ISA ./ T);

v = M.*a;
ddt_h = 0;
ddt_v = 0;
alpha = 1.11*(rho/rho_SL)-0.11;
beta = 1;

T2W_Cruise = Master_Eqn(CD_0,k1,k2,CD_R,n,v,rho,ddt_h,ddt_v,alpha,beta,W2S);

%%%%%%%% Climb %%%%%%%%

alt = 0;
n = 1;
M = 0.8;

[T_ISA,a_ISA,~,rho_ISA] = atmosisa(alt,"extended","on");
T = T_ISA + DeltaT;
a = sqrt(gamma*R*T);
rho = rho_ISA .* (T_ISA ./ T);

v = M.*a;
ddt_h = 28000*0.3048/60;
ddt_v = 0;
alpha = 1.11*(rho/rho_SL)-0.11;
beta = 1.0;

T2W_Climb = Master_Eqn(CD_0,k1,k2,CD_R,n,v,rho,ddt_h,ddt_v,alpha,beta,W2S);

%%%%%%% Ceiling %%%%%%%%

alt = 38700*0.3048;
n = 1;
M = 0.9;

[T_ISA,a_ISA,~,rho_ISA] = atmosisa(alt,"extended","on");
T = T_ISA + DeltaT;
a = sqrt(gamma*R*T);
rho = rho_ISA .* (T_ISA ./ T);

v = M.*a;
ddt_h = 200*0.3048/60;
ddt_v = 0;
alpha = 1.11*(rho/rho_SL)-0.11;
beta = 0.7;

T2W_Ceiling = Master_Eqn(CD_0,k1,k2,CD_R,n,v,rho,ddt_h,ddt_v,alpha,beta,W2S);

%%%%%%% Sustained Turn %%%%%%%%

alt = 20000*0.3048;
n = 4.32;
v = 295;

[T_ISA,a_ISA,~,rho_ISA] = atmosisa(alt,"extended","on");
T = T_ISA + DeltaT;
rho = rho_ISA .* (T_ISA ./ T);

ddt_h = 0;
ddt_v = 0;
alpha = 1.11*(rho/rho_SL)-0.11;
beta = 0.7;

T2W_SustainedTurn = Master_Eqn(CD_0,k1,k2,CD_R,n,v,rho,ddt_h,ddt_v,alpha,beta,W2S);

%%%%%%% Max Speed %%%%%%%%

alt = 30000*0.3048;
n = 1;
M = 1.6;

[T_ISA,a_ISA,~,rho_ISA] = atmosisa(alt,"extended","on");
T = T_ISA + DeltaT;
a = sqrt(gamma*R*T);
rho = rho_ISA .* (T_ISA ./ T);

v = M.*a;
ddt_h = 0;
ddt_v = 0;
alpha = 1.11*(rho/rho_SL)-0.11;
beta = 0.7;

T2W_MaxSpeed = Master_Eqn(CD_0,k1,k2,CD_R,n,v,rho,ddt_h,ddt_v,alpha,beta,W2S);

%% Takeoff Constraint Analysis

alt = 0;
[~,~,rho] = atmos(alt); % already density-scaled
Cl = 1.7/cosd(24);

T2W_Takeoff = TakeoffConstraintAnalysis(TOP,rho/1.225,Cl,W2S);

%% Vertical Constraint Analysis

%%%%%%% Landing %%%%%%%%

alt = 0;
n = 1;
beta = 0.56;
Cl = 1.8/cosd(24);
k = 1.2;
v = 131*1.1*0.514444;

[T_ISA,~,~,rho_ISA] = atmosisa(alt,"extended","on");
T = T_ISA + DeltaT;
rho = rho_ISA .* (T_ISA ./ T);

W2S_Landing = verticalConstraintAnalysis(n,beta,Cl,k,v,rho);

%%%%%%% Stall %%%%%%%%

alt = 0; % Altitude {m}
n = 1; % Load Factor
beta = 0.56; % Weight Lapse
Cl = 1.2; % Lift Coefficient
k = 1.0; % Safety Factor
v = 131*0.514444; % Velocity {m/s}
[T_ISA,~,~,rho_ISA] = atmosisa(alt,"extended","on");
T = T_ISA + DeltaT;
rho = rho_ISA .* (T_ISA ./ T);
W2S_Stall = verticalConstraintAnalysis(n,beta,Cl,k,v,rho);

%% Plotting

figure
hold on
ax = gca;
ax.ColorOrder = lines(13);
set(ax,'DefaultLineLineWidth',1.2)
box on

plot(W2S.*0.02088547,T2W_Cruise)
plot(W2S.*0.02088547,T2W_Climb)
plot(W2S.*0.02088547,T2W_Ceiling)
plot(W2S.*0.02088547,T2W_SustainedTurn)
plot(W2S.*0.02088547,T2W_MaxSpeed)
plot(W2S.*0.02088547,T2W_Takeoff)

xline(W2S_Landing.*0.02088547,'-r')
xline(W2S_Stall.*0.02088547,'-b')

plot(115,0.75,'o')

xlabel("Wing Loading (lb/ft^2)")
ylabel("Thrust-to-Weight")
ylim([0 2.5])

legend("Cruise","Climb","Ceiling","Turn","Max Speed","Takeoff","Landing","Stall","Point")
