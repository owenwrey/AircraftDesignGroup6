% Calculating aircraft CG - takeoff and beta
% Prof. Chakraborty - 
% fall 2025
%--------------------------------------------------------------------------
clc; clear; close all

% data input: locWeight = [X_loc weight]
addpath(genpath('CGFunctions'))

% given geometry

n_lim = 3.8;      % limit load
Wto = 76965;      % MTOW, lbs
tc = .06;      % thickness to chord
b = 48.66; 
bw = b/2;


lam = 1/4;        % wing taper ratio

AR = 4;         % aspect ratio
S = 500;        % wing area

Ymac = 9.73;
Croot = 19.47;
MAC = 13.63;

lam_LE = atand(tand(24) + (0.25*(1+lam)*bw)/((1-lam)*Croot));     % leading edge sweep
%lam_LE = 24; 

% given weights
Wavionic = 2500; % constraint

WJdam = 1015; % wikipedia
Wa120 = 356; %lb - wiki
Wa9 = 188; %lb - wiki

figure
hold on
axis([0 50 0 5])
plot([0 50], [0 0],DisplayName="lengthfuse")
plot(linspace(50-154/12,50,3), [1 1 1],'-o', DisplayName="engine")
plot(linspace(5,6,3), [1 1 1],'-o', DisplayName="avionics1")
plot([6 9 11], [1 1 1],'-o', DisplayName="cockpit")
plot(linspace(11,15,3), [1 1 1],'-o', DisplayName="avionics2")
plot(linspace(15,25,3), [1 1 1],'-o', DisplayName="Fuel1")
plot(linspace(25,35,3), [1 1 1],'-o', DisplayName="Fuel2")
title("Weight Distribution")
xlabel("length along fuselage (ft)")

legend



%Wingless CG

locweights = [];

%avionics 1
locweights(1,:) = [5.5, 500]; % from mission requirements

%pilot
locweights(2,:) = [9, 300]; % from weight conv script

%avionics 2
locweights(3,:) = [13, 2000];

%-fuels--------------------------------------------------------------------
% fuselage diameter = 5.4 ft
% tank dimensions - L = 10ft, ,D = 5 ft
% density jp-5 49.9423696268488 lb/ft^3
volTank = 10*pi*(5/2)^2;
Djp5 = 49.9423696268488;
Wtank = volTank*Djp5;

% fuel tanks
for i = 4:5
locweights(i,:) = [20+10*(i-4),Wtank];
end
%--------------------------------------------------------------------------

% engine
locweights(6,:) = [50-154/24, 2450]; % https://www.mtu.de/engines/military-aircraft-engines/fighter-aircraft/f414/

% fuselage
locweights(7,:) = [50*0.45, Wfuse(Wto,50,5.4)];
plot(50*0.45,0.5,'*',DisplayName="fuselage")

%CG
CGnowing = (sum(locweights(:,1).*locweights(:,2)))./sum(locweights(:,2));
plot(CGnowing,1.5,'*',DisplayName="CGnowing")


% main wing
% 30% of MAC, 
wingG = polyshape([0 19.248 (48.66/2)*tand(lam_LE)+ 4.8675 (48.66/2)*tand(lam_LE)],[0 0 (48.66/2) (48.66/2)]);
[Gwing,ygw] = centroid(wingG);

Ww = Wwing(n_lim,Wto,tc,lam_LE,lam,AR,S);
trpMAC =(9.73)*tand(lam_LE)+0.3*MAC;
Xw = CGnowing-trpMAC+Gwing;

plot(Xw,2,'*',DisplayName="cgwing")
plot(CGnowing,2,'*',DisplayName="30%MAC")
plot([CGnowing-trpMAC CGnowing-trpMAC+Croot],[2 2],DisplayName="wingroot")

% plots wing
if false
figure
plot(wingG)
hold on
plot(Gwing,ygw,"*")
xline(trpMAC)
yline(9.73)
axis([0 30 0 30])
end

locweights(end+1,:) = [Xw, Ww]; %add to weight buildup

% horizontal tails
%lam_LEH = atand(tand(30) + (0.25*(1+0.4)*(22.17/2))/((1-0.4)*7.92)); 
lam_LEH = 30;

HstabG = polyshape([0 7.92 (22.17/2)*tand(lam_LEH)+ 3.17 (22.17/2)*tand(lam_LEH)],[0 0 (22.17/2) (22.17/2)]);
[Ghstab,~] = centroid(HstabG);

WHs = Whstab(Wto,n_lim,122.93,22.17/2,.231,5.88,0.525*50);

xqchH = 0.525*50 + CGnowing-trpMAC+0.25*MAC;

XHs = xqchH - 0.25*5.88-(4.75)*tand(lam_LEH) + Ghstab;
plot([XHs-Ghstab,XHs-Ghstab+7.92], [2,2],DisplayName="hStabRoot")
plot(xqchH,2,"*",DisplayName="H qchord")
plot(XHs,2,"*",DisplayName="H C.G.")



%plot([],[2 2],'*',DisplayName="cgwing")

if false
figure
hold on
plot(HstabG)
axis([0 20 0 20])
xline(4.75*tand(lam_LEH) + 0.25*5.88)
yline(4.75)
end

locweights(end+1,:) = [XHs, WHs];

% Vert Stab
% Ybar = 3.04


lam_LEV = atand(tand(35) + (0.25*(1+lam)*bw)/((1-lam)*Croot));
%lam_LEV = 35;

VstabG = polyshape([0 11.37 (7.39)*tand(lam_LEV)+ 3.41 (7.39)*tand(lam_LEV)],[0 0 (7.39) (7.39)]);
[Gvstab,yvs] = centroid(VstabG);

WVs = 2*Wvstab(7.39,Wto,n_lim,54.59,0.2*54.59,1.7,.475*50,1,.3,35);

xqcVs = 0.475*50 + CGnowing-trpMAC+0.25*MAC;

XVs = xqcVs - 0.25*8.10-3.04*tand(lam_LEV) + Gvstab;
plot([XVs-Gvstab,XVs-Gvstab+11.37], [2.5,2.5],DisplayName="vStabRoot")
plot(xqcVs,2.5,"*",DisplayName="V qchord")
plot(XVs,2.5,"*",DisplayName="V C.G.")

if false
figure
hold on
plot(VstabG)
plot(Gvstab,yvs,"*")
axis([0 12 0 12])
yline(2*1.52)
end

locweights(end+1,:) = [XVs, WVs];

% jdams + aim9s

bbomb = linspace(0,b/2,4);
bbomb = bbomb(2:3);

cbomb = linspace(Croot/2,(28.36-23.5)/2+23.5,4);
cbomb = cbomb(2:3);

locweights([end+1 end+2],:) = [CGnowing-trpMAC+cbomb', 2*WJdam*[1;1]];

locweights(end+1,:) = [CGnowing-trpMAC+(28.36-23.5)/2, 2*Wa9];


%-Add LG math here---------------------------------------------------------






%--------------------------------------------------------------------------

% Final CG

CG = (sum(locweights(:,1).*locweights(:,2)))./sum(locweights(:,2));

plot(CG,3,".",DisplayName="Center of G.",MarkerSize=24)
xline(CG)

CG

locweights
