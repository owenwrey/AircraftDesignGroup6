clc
% clear
% segment names
SegNames = {'SWT','TKO','CLIMB','CR OBD', 'DESC 1', 'COMBAT', 'WP FIRE', 'CR IBD', ...
'DESC 2', 'LTR2', 'LTS'};
% points per segment
npts = [...
3; % 1. start, warmup, taxi (SWT)
3; % 2. takeoff
10; % 3. climb 1
10; % 4. cruise outbound
10; % 5. descent 1
5; % 6. combat
1; % 7. drop/fire ordnance
10; % 8. climb 2
10; % 9. cruise inbound
10; % 10. descend 2
10; % 11. loiter
3; % 12. landing, taxi, shutdown (LTS)
];
npts_sum = sum(npts);
Tbl = table();
Tbl.SegNum = zeros(npts_sum,1);
Tbl.Seg = repmat({' '},npts_sum,1);
ctr = 1;
for i = 1:length(SegNames)
Tbl.Seg([ctr:ctr+npts(i)-1]) = SegNames(i);
Tbl.SegNum([ctr:ctr+npts(i)-1]) = i;
ctr = ctr + npts(i);
end
Tbl.Time = zeros(npts_sum,1); % time (min)
Tbl.dTime = zeros(npts_sum,1); % delta time (min)
Tbl.Alt = zeros(npts_sum,1); % altitude (ft)
Tbl.rho = zeros(npts_sum,1); % density (slug/ft^3)
Tbl.KEAS = zeros(npts_sum,1); % equivalent airspeed (kt)
Tbl.KTAS = zeros(npts_sum,1); % true airspeed (kt)
Tbl.MACH = zeros(npts_sum,1); % Mach number
Tbl.GS = zeros(npts_sum,1); % ground speed (kt)
Tbl.FPA = zeros(npts_sum,1); % flightpath angle (deg)
Tbl.EnHt = zeros(npts_sum,1); % energy height (ft)
Tbl.Dist = zeros(npts_sum,1); % distance (NM)
Tbl.dDist = zeros(npts_sum,1); % delta distance (NM)
Tbl.dhdt = zeros(npts_sum,1); % rate of climb (ft/min)
Tbl.dVdt = zeros(npts_sum,1); % acceleration
Tbl.Weight = zeros(npts_sum,1); % weight (lb)
Tbl.WtDrop = zeros(npts_sum,1); % dropped weight (lb)
Tbl.WtFrac = zeros(npts_sum,1); % weight fraction
Tbl.CL = zeros(npts_sum,1); % lift coefficient
Tbl.CD0 = zeros(npts_sum,1); % drag polar
Tbl.K1 = zeros(npts_sum,1); % drag polar
Tbl.K2 = zeros(npts_sum,1); % drag polar
Tbl.CDR = zeros(npts_sum,1); % drag polar
Tbl.CD = zeros(npts_sum,1); % drag coefficient
Tbl.L_D = zeros(npts_sum,1); % lift-to-drag ratio
Tbl.Drag = zeros(npts_sum,1); % drag (lbf)
Tbl.Thrust = zeros(npts_sum,1); % thrust (lbf)
Tbl.TLapse = zeros(npts_sum,1); % thrust lapse
Tbl.Ps = zeros(npts_sum,1); % specific excess pwr (ft/min)
Tbl.THROT = zeros(npts_sum,1); % throttle setting
Tbl.FF = zeros(npts_sum,1); % fuel flow (lb/h)
Tbl.dFuel = zeros(npts_sum,1); % delta fuel (lb)
Tbl.FuelRem = zeros(npts_sum,1); % fuel remaining (lb)
Tbl.FuelBurn = zeros(npts_sum,1); % fuel burned (lb)


TLapse = griddedInterpolant([0, 10000, 20000, 30000, 40000, 50000], [1, 0.80, 0.60, ...
0.40, 0.20, 0.15],'linear','nearest');
Alt_ft = 28000;
TLapse(Alt_ft);


W0 = 63738;
FuelReq = 16440;
tol = 2;

while tol >1

W_S = 130; % takeoff wing loading (psf)
rho_SL = 0.0023769;
mps2kts = 1.94384; % meters per sec to knots
kts2fps = 1/0.59248; % knots to feet per sec
NM2ft = 6067; % nautical miles to feet
ft2m = 0.305; % feet to meters
Thrust = 44000;
% segment information

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. start, warmup, taxi
% 3 min duration
% no distance credit

% given

for i = 2:npts(1,:)
SFC_idle = 0.6; %lb/(lb.h)
T_W_idle = 0.05;


Tbl.Time(1:npts(1,:)) = linspace(0,3,npts(1,:)); % time (min)
Tbl.dTime(1) = 0;
Tbl.dTime(i) = Tbl.Time(i) - Tbl.Time(i-1); % delta time (min)
Tbl.Alt(i) = 0; % altitude (ft)
[T, a, P, rho] = atmosisa(Tbl.Alt(i)*ft2m);
Tbl.rho(1) = rho*0.00194032;
Tbl.rho(i) = rho*0.00194032; % density in slug/ft^3
Tbl.KEAS(i) = 0; % equivalent airspeed (kt)
Tbl.KTAS(i) = 0; % true airspeed (kt)
Tbl.MACH(i) = 0; % Mach number
Tbl.GS(i) = 0; % ground speed (kt)
Tbl.FPA(i) = 0; % flightpath angle (deg)
Tbl.Dist(i) = 0; % distance (NM)1
Tbl.dDist(i) = 0; % delta distance (NM)
Tbl.dhdt(i) = 0; % rate of climb (ft/min)
Tbl.dVdt(i) = 0; % acceleration
Tbl.CL(i) = 0; % lift coefficient
Tbl.CD0(i) = 0.008; % drag polar
Tbl.K1(i) = 0; % drag polar
Tbl.K2(i) = 0.08; % drag polar
Tbl.CDR(i) = 0; % drag polar
Tbl.CD(i) = 0; % drag coefficient
Tbl.L_D(i) = 0; % lift-to-drag ratio
Tbl.Drag(i) = 0; % drag (lbf)
Tbl.Thrust(1) = T_W_idle*W0;
Tbl.Thrust(i) = T_W_idle*W0; % thrust (lbf)
Tbl.TLapse(1) = 1;
Tbl.TLapse(i) = 1; % thrust lapse
Tbl.Ps(i) = 0; % specific excess pwr (ft/min)
Tbl.THROT(i) = 0; % throttle setting
Tbl.FF(1) = SFC_idle*Tbl.Thrust(1);
Tbl.FF(i) = SFC_idle*Tbl.Thrust(i); % fuel flow (lb/h)
Tbl.dFuel(i) = (Tbl.FF(i)*Tbl.dTime(i))/60; % delta fuel (lb)
Tbl.FuelBurn(1) = Tbl.dFuel(1);
Tbl.FuelBurn(i) = sum(Tbl.dFuel(1:i)); % fuel burned (lb)
Tbl.FuelRem(1) = FuelReq - Tbl.dFuel(1);
Tbl.FuelRem(i) = Tbl.FuelRem(i-1) - Tbl.dFuel(i); % fuel remaining (lb)
Tbl.Weight(1) = W0;
Tbl.Weight(i) = Tbl.Weight(i-1) - Tbl.dFuel(i); % weight (lb)
Tbl.WtDrop(i) = 0; % dropped weight (lb)
Tbl.WtFrac(1) = 1;
Tbl.WtFrac(i) = Tbl.Weight(i)/Tbl.Weight(i-1);% weight fraction
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. takeoff
% 1 min duration
% full afterburner (THROT = 1.25)
% no distance credit
istart = npts(1,:)+1;
iend = sum(npts(1:2));

SFC_takeoff = 1.85;

for  i = istart:iend

Tbl.Time(istart:iend) = linspace(3,4,npts(2,:)); % time (min)
Tbl.dTime(i) = Tbl.Time(i) - Tbl.Time(i-1);% delta time (min)
Tbl.Alt(i) = 0; % altitude (ft)
[T, a, P, rho] = atmosisa(Tbl.Alt(i)*ft2m);
Tbl.rho(i) = rho*0.00194032; % density in slug/ft^3
Tbl.KEAS(i) = 0; % equivalent airspeed (kt)
Tbl.KTAS(i) = 0; % true airspeed (kt)
Tbl.MACH(i) = 0; % Mach number
Tbl.GS(i) = 0; % ground speed (kt)
Tbl.FPA(i) = 0; % flightpath angle (deg)
Tbl.Dist(i) = 0; % distance (NM)
Tbl.dDist(i) = 0; % delta distance (NM)
Tbl.dhdt(i) = 0; % rate of climb (ft/min)
Tbl.dVdt(i) = 0; % acceleration
Tbl.CL(i) = 0; % lift coefficient
Tbl.CD0(i) = 0.008; % drag polar
Tbl.K1(i) = 0; % drag polar
Tbl.K2(i) = 0.08; % drag polar
Tbl.CDR(i) = 0; % drag polar
Tbl.CD(i) = 0; % drag coefficient
Tbl.L_D(i) = 0; % lift-to-drag ratio
Tbl.Drag(i) = 0; % drag (lbf)
Tbl.TLapse(i) = 1; % thrust lapse
Tbl.Ps(i) = 0; % specific excess pwr (ft/min)
Tbl.THROT(i) = 1.25; % throttle setting
Tbl.Thrust(i) = Thrust*Tbl.TLapse(i)*Tbl.THROT(i); % thrust (lbf)
Tbl.FF(i) = SFC_takeoff*Tbl.Thrust(i); % fuel flow (lb/h)
Tbl.dFuel(i) = (Tbl.FF(i)*Tbl.dTime(i))/60; % delta fuel (lb)
Tbl.FuelBurn(i) = sum(Tbl.dFuel(1:i)); % fuel burned (lb)
Tbl.FuelRem(i) = Tbl.FuelRem(i-1) - Tbl.dFuel(i); % fuel remaining (lb)
Tbl.Weight(i) = Tbl.Weight(i-1) - Tbl.dFuel(i); % weight (lb)
Tbl.WtDrop(i) = 0; % dropped weight (lb)
Tbl.WtFrac(i) = Tbl.Weight(i)/Tbl.Weight(i-1);% weight fraction
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 3. climb 1
% Mach 0.8
% level off at 25,000 ft
% no afterburner [THROT = 1]
% take distance credit

istart3 = sum(npts(1:2))+1;
iend3 = sum(npts(1:3));
SFC_climb = 0.75;
for  i = istart3:iend3

Tbl.Alt(istart3:iend3) = linspace(0,25000,npts(3,:)); % altitude (ft)

[T, a, P, rho] = atmosisa(Tbl.Alt(i)*ft2m);
Tbl.rho(i) = rho*0.00194032; % density in slug/ft^3
Tbl.TLapse(i) = TLapse(Tbl.Alt(i)); % thrust lapse
Tbl.MACH(i) = 0.8; % Mach number
Tbl.KTAS(i) = Tbl.MACH(i)*a*mps2kts; % true airspeed (kt)
Tbl.KEAS(i) = Tbl.KTAS(i)*sqrt(Tbl.rho(i)/rho_SL); % equivalent airspeed (kt)

Tbl.THROT(i) = 1; % throttle setting
Tbl.Thrust(i) = Thrust*Tbl.TLapse(i)*Tbl.THROT(i); % thrust (lbf)
Tbl.FF(i) = SFC_climb*Tbl.Thrust(i); % fuel flow (lb/h)
Tbl.dhdt(i) = ((Tbl.Thrust(i)*Tbl.KTAS(i)*kts2fps*60) - (Tbl.Drag(i)*Tbl.KTAS(i)*kts2fps*60))/Tbl.Weight(i-1) ; % rate of climb (ft/min)
Tbl.Ps(i) = Tbl.dhdt(i); % specific excess pwr (ft/min)

Tbl.CL(i) = (Tbl.WtFrac(i-1)*W_S)/(0.5*Tbl.rho(i)*(Tbl.KTAS(i)*kts2fps)^2); % lift coefficient
Tbl.CD0(i) = 0.008; % drag polar
Tbl.K1(i) = 0; % drag polar
Tbl.K2(i) = 0.08; % drag polar
Tbl.CDR(i) = 0; % drag polar
Tbl.CD(i) = Tbl.CD0(i) + Tbl.K1(i)*Tbl.CL(i) + Tbl.K2(i)*Tbl.CL(i)^2 + Tbl.CDR(i); % drag coefficient
Tbl.L_D(i) = Tbl.CL(i)/Tbl.CD(i); % lift-to-drag ratio
Tbl.Drag(i) = Tbl.CD(i) * 0.5* Tbl.rho(i) * (Tbl.KTAS(i)*kts2fps)^2 * (W0/W_S); % drag (lbf)

Tbl.dTime(i) = (Tbl.Alt(i) - Tbl.Alt(i-1))/Tbl.dhdt(i);% delta time (min)
Tbl.Time(i) = Tbl.dTime(i)+Tbl.Time(i-1); % time (min)
Tbl.GS(i) = Tbl.KTAS(i)*sind(Tbl.FPA(i)); % ground speed (kt)
Tbl.FPA(i) = asind(Tbl.dhdt(i)/(Tbl.KTAS(i)*kts2fps*60)); % flightpath angle (deg)
Tbl.Dist(i) = (Tbl.Alt(i) - Tbl.Alt(i-1))/(tand(Tbl.FPA(i))*NM2ft); % distance (NM)
Tbl.dDist(i) =Tbl.Dist(i) - Tbl.Dist(i-1); % delta distance (NM)
Tbl.dVdt(istart3) = 0;
Tbl.dVdt(i) = (Tbl.KTAS(i)*kts2fps - (Tbl.KTAS(i-1)*kts2fps))/(Tbl.dTime(i)*60); % acceleration


Tbl.dFuel(i) = (Tbl.FF(i)*Tbl.dTime(i))/60; % delta fuel (lb)
Tbl.FuelBurn(i) = sum(Tbl.dFuel(1:i)); % fuel burned (lb)
Tbl.FuelRem(i) = Tbl.FuelRem(i-1) - Tbl.dFuel(i); % fuel remaining (lb)
Tbl.Weight(i) = Tbl.Weight(i-1) - Tbl.dFuel(i); % weight (lb)
Tbl.WtDrop(i) = 0; % dropped weight (lb)
Tbl.WtFrac(i) = Tbl.Weight(i)/Tbl.Weight(i-1);% weight fraction

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 4. cruise outbound

% 25,000 ft
% 750 NM
SFC_cruise = 0.85; %estimate
istart4 = sum(npts(1:3))+1;
iend4 = sum(npts(1:4));

for  i = istart4:iend4
Tbl.Dist(istart4:iend4) = linspace((Tbl.Dist(istart4-1)),(Tbl.Dist(istart4-1)+ 750),npts(4,:)); % distance (NM)
Tbl.dDist(i) = Tbl.Dist(i) - Tbl.Dist(i-1); % delta distance (NM)
Tbl.Alt(i) = 25000; % Altitude (ft)
[T, a, P, rho] = atmosisa(Tbl.Alt(i)*ft2m);
Tbl.rho(i) = rho*0.00194032; % density in slug/ft^3
Tbl.TLapse(i) = TLapse(Tbl.Alt(i)); % thrust lapse
Tbl.KTAS(i) = 662; % true airspeed (kt) (660 mph taken from F18 e/f cruise)
Tbl.KEAS(i) = Tbl.KTAS(i)*sqrt(Tbl.rho(i)/rho_SL); % equivalent airspeed (kt)
Tbl.MACH(i) = Tbl.KTAS(i)/(a*mps2kts); % Mach number

Tbl.CL(i) = (Tbl.WtFrac(i-1)*W_S)/(0.5*Tbl.rho(i)*(Tbl.KTAS(i)*kts2fps)^2); % lift coefficient
Tbl.CD0(i) = 0.008; % drag polar
Tbl.K1(i) = 0; % drag polar
Tbl.K2(i) = 0.08; % drag polar
Tbl.CDR(i) = 0.01; % drag polar
Tbl.CD(i) = Tbl.CD0(i) + Tbl.K1(i)*Tbl.CL(i) + Tbl.K2(i)*Tbl.CL(i)^2 + Tbl.CDR(i); % drag coefficient
Tbl.L_D(i) = Tbl.CL(i)/Tbl.CD(i); % lift-to-drag ratio
Tbl.Drag(i) = Tbl.CD(i) * 0.5* Tbl.rho(i) * (Tbl.KTAS(i)*kts2fps)^2 * (W0/W_S); % drag (lbf)

% solve for throttle setting 
Tbl.Thrust(i) = Tbl.Drag(i); % thrust (lbf)
Tbl.THROT(i) = Tbl.Thrust(i)/(Thrust*Tbl.TLapse(i)); % throttle setting
Tbl.FF(i) = SFC_cruise*Tbl.Thrust(i); % fuel flow (lb/h)
Tbl.dhdt(i) = 0 ; % rate of climb (ft/min)
Tbl.Ps(i) = 0; % specific excess pwr (ft/min)

Tbl.dTime(i) = Tbl.dDist(i)/(Tbl.KTAS(i)/60);% delta time (min)
Tbl.Time(i) = Tbl.dTime(i)+Tbl.Time(i-1); % time (min)
Tbl.FPA(i) = 0; % flightpath angle (deg)
Tbl.GS(i) = Tbl.KTAS(i)*sind(Tbl.FPA(i)); % ground speed (kt)
Tbl.dVdt(i) = 0; % acceleration

Tbl.dFuel(istart4) = 0;
Tbl.dFuel(i) = (Tbl.FF(i)*Tbl.dTime(i))/60; % delta fuel (lb)
Tbl.FuelBurn(i) = sum(Tbl.dFuel(1:i)); % fuel burned (lb)
Tbl.FuelRem(i) = Tbl.FuelRem(i-1) - Tbl.dFuel(i); % fuel remaining (lb)
Tbl.Weight(i) = Tbl.Weight(i-1) - Tbl.dFuel(i); % weight (lb)
Tbl.WtDrop(i) = 0; % dropped weight (lb)
Tbl.WtFrac(i) = Tbl.Weight(i)/Tbl.Weight(i-1);% weight fraction
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 5. descent 1
% -3,000 ft/min

% finish at sea level
SFC_descent = 0.60;
istart5 = sum(npts(1:4))+1;
iend5 = sum(npts(1:5));
for  i = istart5:iend5

Tbl.Alt(istart5:iend5) = linspace(Tbl.Alt(istart5 - 1),0,npts(5,:)); % altitude (ft)

[T, a, P, rho] = atmosisa(Tbl.Alt(i)*ft2m);
Tbl.rho(i) = rho*0.00194032; % density in slug/ft^3
Tbl.TLapse(i) = TLapse(Tbl.Alt(i)); % thrust lapse
Tbl.MACH(i) = 0.8; % Mach number
Tbl.KTAS(i) = Tbl.MACH(i)*a*mps2kts; % true airspeed (kt)
Tbl.KEAS(i) = Tbl.KTAS(i)*sqrt(Tbl.rho(i)/rho_SL); % equivalent airspeed (kt)

Tbl.THROT(i) = 1; % throttle setting
Tbl.Thrust(i) = Thrust*Tbl.TLapse(i)*Tbl.THROT(i); % thrust (lbf)
Tbl.FF(i) = SFC_descent*Tbl.Thrust(i); % fuel flow (lb/h)
Tbl.dhdt(i) = -3000; % rate of climb (ft/min)
Tbl.Ps(i) = Tbl.dhdt(i); % specific excess pwr (ft/min)

Tbl.CL(i) = (Tbl.WtFrac(i-1)*W_S)/(0.5*Tbl.rho(i)*(Tbl.KTAS(i)*kts2fps)^2); % lift coefficient
Tbl.CD0(i) = 0.008; % drag polar
Tbl.K1(i) = 0; % drag polar
Tbl.K2(i) = 0.08; % drag polar
Tbl.CDR(i) = 0; % drag polar
Tbl.CD(i) = Tbl.CD0(i) + Tbl.K1(i)*Tbl.CL(i) + Tbl.K2(i)*Tbl.CL(i)^2 + Tbl.CDR(i); % drag coefficient
Tbl.L_D(i) = Tbl.CL(i)/Tbl.CD(i); % lift-to-drag ratio
Tbl.Drag(i) = Tbl.CD(i) * 0.5* Tbl.rho(i) * (Tbl.KTAS(i)*kts2fps)^2 * (W0/W_S); % drag (lbf)

Tbl.dTime(i) = (Tbl.Alt(i) - Tbl.Alt(i-1))/Tbl.dhdt(i);% delta time (min)
Tbl.Time(i) = Tbl.dTime(i)+Tbl.Time(i-1); % time (min)
Tbl.GS(i) = Tbl.KTAS(i)*sind(Tbl.FPA(i)); % ground speed (kt)
Tbl.FPA(i) = asind(Tbl.dhdt(i)/(Tbl.KTAS(i)*kts2fps*60)); % flightpath angle (deg)
Tbl.Dist(istart5) = Tbl.Dist(istart5 - 1);
Tbl.Dist(i) = Tbl.GS(i)*Tbl.dTime(i)/60 + Tbl.Dist(istart5 - 1); % distance (NM)
Tbl.dDist(i) =Tbl.Dist(i) - Tbl.Dist(i-1); % delta distance (NM)
Tbl.dVdt(i) = (Tbl.KTAS(i)*kts2fps - (Tbl.KTAS(i-1)*kts2fps))/(Tbl.dTime(i)*60); % acceleration
Tbl.dVdt(istart5) = 0;

Tbl.dFuel(i) = (Tbl.FF(i)*Tbl.dTime(i))/60; % delta fuel (lb)
Tbl.FuelBurn(i) = sum(Tbl.dFuel(1:i)); % fuel burned (lb)
Tbl.FuelRem(i) = Tbl.FuelRem(i-1) - Tbl.dFuel(i); % fuel remaining (lb)
Tbl.Weight(i) = Tbl.Weight(i-1) - Tbl.dFuel(i); % weight (lb)
Tbl.WtDrop(i) = 0; % dropped weight (lb)
Tbl.WtFrac(i) = Tbl.Weight(i)/Tbl.Weight(i-1);% weight fraction

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 6. combat
% take distance credit - 100 nm
% full afterburner
% sea level
% Mach 0.85

SFC_combat = 0.85; % reevaluate
istart6 = sum(npts(1:5))+1;
iend6 = sum(npts(1:6));

for  i = istart6:iend6
Tbl.Dist(istart6:iend6) = linspace((Tbl.Dist(istart6-1)),(Tbl.Dist(istart6-1)+ 100),npts(6,:)); % distance (NM)
Tbl.dDist(i) = Tbl.Dist(i) - Tbl.Dist(i-1); % delta distance (NM)
Tbl.Alt(i) = 0; % Altitude (ft)
[T, a, P, rho] = atmosisa(Tbl.Alt(i)*ft2m);
Tbl.rho(i) = rho*0.00194032; % density in slug/ft^3
Tbl.TLapse(i) = TLapse(Tbl.Alt(i)); % thrust lapse
Tbl.KTAS(i) = 562; % true airspeed (kt) (660 mph taken from F18 e/f cruise)
Tbl.KEAS(i) = Tbl.KTAS(i)*sqrt(Tbl.rho(i)/rho_SL); % equivalent airspeed (kt)
Tbl.MACH(i) = Tbl.KTAS(i)/(a*mps2kts); % Mach number

Tbl.CL(i) = (Tbl.WtFrac(i-1)*W_S)/(0.5*Tbl.rho(i)*(Tbl.KTAS(i)*kts2fps)^2); % lift coefficient
Tbl.CD0(i) = 0.008; % drag polar
Tbl.K1(i) = 0; % drag polar
Tbl.K2(i) = 0.08; % drag polar
Tbl.CDR(i) = 0.01; % drag polar
Tbl.CD(i) = Tbl.CD0(i) + Tbl.K1(i)*Tbl.CL(i) + Tbl.K2(i)*Tbl.CL(i)^2 + Tbl.CDR(i); % drag coefficient
Tbl.L_D(i) = Tbl.CL(i)/Tbl.CD(i); % lift-to-drag ratio
Tbl.Drag(i) = Tbl.CD(i) * 0.5* Tbl.rho(i) * (Tbl.KTAS(i)*kts2fps)^2 * (W0/W_S); % drag (lbf)

% solve for throttle setting 
Tbl.Thrust(i) = Tbl.Drag(i); % thrust (lbf)
Tbl.THROT(i) = Tbl.Thrust(i)/(Thrust*Tbl.TLapse(i)); % throttle setting
Tbl.FF(i) = SFC_cruise*Tbl.Thrust(i); % fuel flow (lb/h)
Tbl.dhdt(i) = 0 ; % rate of climb (ft/min)
Tbl.Ps(i) = 0; % specific excess pwr (ft/min)

Tbl.dTime(i) = Tbl.dDist(i)/(Tbl.KTAS(i)/60);% delta time (min)
Tbl.Time(i) = Tbl.dTime(i)+Tbl.Time(i-1); % time (min)
Tbl.FPA(i) = 0; % flightpath angle (deg)
Tbl.GS(i) = Tbl.KTAS(i)*sind(Tbl.FPA(i)); % ground speed (kt)
Tbl.dVdt(i) = 0; % acceleration

Tbl.dFuel(istart) = 0;
Tbl.dFuel(i) = (Tbl.FF(i)*Tbl.dTime(i))/60; % delta fuel (lb)
Tbl.FuelBurn(i) = sum(Tbl.dFuel(1:i)); % fuel burned (lb)
Tbl.FuelRem(i) = Tbl.FuelRem(i-1) - Tbl.dFuel(i); % fuel remaining (lb)
Tbl.Weight(i) = Tbl.Weight(i-1) - Tbl.dFuel(i); % weight (lb)
Tbl.WtDrop(i) = 0; % dropped weight (lb)
Tbl.WtFrac(i) = Tbl.Weight(i)/Tbl.Weight(i-1);% weight fraction
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 7. weapons fire/drop
%  4380 lb dropped
% occurs only in M1
SFC_weightdrop = 0.85;
istart7 = sum(npts(1:6))+1;
iend7 = sum(npts(1:7));

for  i = istart7:iend7
Tbl.Time(i) = Tbl.Time(istart7 - 1); % distance (NM)
Tbl.dTime(i) = 0; % delta distance (NM)
Tbl.Dist(i) = Tbl.Dist(istart7 - 1);
Tbl.dDist(i) = 0;

Tbl.Alt(i) = 0; % Altitude (ft)
[T, a, P, rho] = atmosisa(Tbl.Alt(i)*ft2m);
Tbl.rho(i) = rho*0.00194032; % density in slug/ft^3
Tbl.TLapse(i) = TLapse(Tbl.Alt(i)); % thrust lapse
Tbl.KTAS(i) = 722.4; % true airspeed (kt)
Tbl.KEAS(i) = Tbl.KTAS(i)*sqrt(Tbl.rho(i)/rho_SL); % equivalent airspeed (kt)
Tbl.MACH(i) = Tbl.KTAS(i)/(a*mps2kts); % Mach number

Tbl.CL(i) = (Tbl.WtFrac(i-1)*W_S)/(0.5*Tbl.rho(i)*(Tbl.KTAS(i)*kts2fps)^2); % lift coefficient
Tbl.CD0(i) = 0.008; % drag polar
Tbl.K1(i) = 0; % drag polar
Tbl.K2(i) = 0.08; % drag polar
Tbl.CDR(i) = 0; % drag polar
Tbl.CD(i) = Tbl.CD0(i) + Tbl.K1(i)*Tbl.CL(i) + Tbl.K2(i)*Tbl.CL(i)^2 + Tbl.CDR(i); % drag coefficient
Tbl.L_D(i) = Tbl.CL(i)/Tbl.CD(i); % lift-to-drag ratio
Tbl.Drag(i) = Tbl.CD(i) * 0.5* Tbl.rho(i) * (Tbl.KTAS(i)*kts2fps)^2 * (W0/W_S); % drag (lbf)

% solve for throttle setting 
Tbl.THROT(i) = 1.25; % throttle setting
Tbl.Thrust(i) = Thrust*Tbl.TLapse(i)*Tbl.THROT(i); % thrust (lbf); % thrust (lbf)
Tbl.FF(i) = SFC_weightdrop*Tbl.Thrust(i); % fuel flow (lb/h)
Tbl.dhdt(i) = 0;  % rate of climb (ft/min)
Tbl.Ps(i) = Tbl.dhdt(i); % specific excess pwr (ft/min)

Tbl.FPA(i) = 0; % flightpath angle (deg)
Tbl.GS(i) = Tbl.KTAS(i)*sind(Tbl.FPA(i)); % ground speed (kt)
Tbl.dVdt(i) = 0; % acceleration

Tbl.dFuel(istart) = 0;
Tbl.dFuel(i) = (Tbl.FF(i)*Tbl.dTime(i))/60; % delta fuel (lb)
Tbl.FuelBurn(i) = sum(Tbl.dFuel(1:i)); % fuel burned (lb)
Tbl.FuelRem(i) = Tbl.FuelRem(i-1) - Tbl.dFuel(i); % fuel remaining (lb)
Tbl.WtDrop(i) = 4380; % dropped weight (lb)
Tbl.Weight(i) = Tbl.Weight(i-1) - Tbl.dFuel(i) - Tbl.WtDrop(i); % weight (lb)
Tbl.WtFrac(i) = 1;% weight fraction
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 8. climb 2
% Mach 0.8
% level off at 25,000 ft
% no afterburner [THROT = 1]
% take distance credit

istart8 = sum(npts(1:7))+1;
iend8 = sum(npts(1:8));
SFC_climb = 0.75;
for  i = istart8:iend8

Tbl.Alt(istart8:iend8) = linspace(0,25000,npts(8,:)); % altitude (ft)

[T, a, P, rho] = atmosisa(Tbl.Alt(i)*ft2m);
Tbl.rho(i) = rho*0.00194032; % density in slug/ft^3
Tbl.TLapse(i) = TLapse(Tbl.Alt(i)); % thrust lapse
Tbl.MACH(i) = 0.8; % Mach number
Tbl.KTAS(i) = Tbl.MACH(i)*a*mps2kts; % true airspeed (kt)
Tbl.KEAS(i) = Tbl.KTAS(i)*sqrt(Tbl.rho(i)/rho_SL); % equivalent airspeed (kt)

Tbl.THROT(i) = 1; % throttle setting
Tbl.Thrust(i) = Thrust*Tbl.TLapse(i)*Tbl.THROT(i); % thrust (lbf)
Tbl.FF(i) = SFC_climb*Tbl.Thrust(i); % fuel flow (lb/h)
Tbl.dhdt(i) = ((Tbl.Thrust(i)*Tbl.KTAS(i)*kts2fps*60) - (Tbl.Drag(i)*Tbl.KTAS(i)*kts2fps*60))/Tbl.Weight(i-1) ; % rate of climb (ft/min)
Tbl.Ps(i) = Tbl.dhdt(i); % specific excess pwr (ft/min)

Tbl.CL(i) = (Tbl.WtFrac(i-1)*W_S)/(0.5*Tbl.rho(i)*(Tbl.KTAS(i)*kts2fps)^2); % lift coefficient
Tbl.CD0(i) = 0.008; % drag polar
Tbl.K1(i) = 0; % drag polar
Tbl.K2(i) = 0.08; % drag polar
Tbl.CDR(i) = 0; % drag polar
Tbl.CD(i) = Tbl.CD0(i) + Tbl.K1(i)*Tbl.CL(i) + Tbl.K2(i)*Tbl.CL(i)^2 + Tbl.CDR(i); % drag coefficient
Tbl.L_D(i) = Tbl.CL(i)/Tbl.CD(i); % lift-to-drag ratio
Tbl.Drag(i) = Tbl.CD(i) * 0.5* Tbl.rho(i) * (Tbl.KTAS(i)*kts2fps)^2 * (W0/W_S); % drag (lbf)

Tbl.dTime(i) = (Tbl.Alt(i) - Tbl.Alt(i-1))/Tbl.dhdt(i);% delta time (min)
Tbl.Time(i) = Tbl.dTime(i)+Tbl.Time(i-1); % time (min)
Tbl.GS(i) = Tbl.KTAS(i)*sind(Tbl.FPA(i)); % ground speed (kt)
Tbl.FPA(i) = asind(Tbl.dhdt(i)/(Tbl.KTAS(i)*kts2fps*60)); % flightpath angle (deg)
Tbl.Dist(i) = (Tbl.Alt(i) - Tbl.Alt(i-1))/(tand(Tbl.FPA(i))*NM2ft); % distance (NM)
Tbl.dDist(i) =Tbl.Dist(i) - Tbl.Dist(i-1); % delta distance (NM)
Tbl.dVdt(istart3) = 0;
Tbl.dVdt(i) = (Tbl.KTAS(i)*kts2fps - (Tbl.KTAS(i-1)*kts2fps))/(Tbl.dTime(i)*60); % acceleration


Tbl.dFuel(i) = (Tbl.FF(i)*Tbl.dTime(i))/60; % delta fuel (lb)
Tbl.FuelBurn(i) = sum(Tbl.dFuel(1:i)); % fuel burned (lb)
Tbl.FuelRem(i) = Tbl.FuelRem(i-1) - Tbl.dFuel(i); % fuel remaining (lb)
Tbl.Weight(i) = Tbl.Weight(i-1) - Tbl.dFuel(i); % weight (lb)
Tbl.WtDrop(i) = 0; % dropped weight (lb)
Tbl.WtFrac(i) = Tbl.Weight(i)/Tbl.Weight(i-1);% weight fraction

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 9. cruise inbound
% 750 nm

% sea level

istart9 = sum(npts(1:8))+1;
iend9 = sum(npts(1:9));

SFC_cruise = 0.85; %estimate


for  i = istart9:iend9
Tbl.Dist(istart9:iend9) = linspace((Tbl.Dist(istart9-1)),(Tbl.Dist(istart9-1)+ 750),npts(9,:)); % distance (NM)
Tbl.dDist(i) = Tbl.Dist(i) - Tbl.Dist(i-1); % delta distance (NM)
Tbl.Alt(i) = 25000; % Altitude (ft)
[T, a, P, rho] = atmosisa(Tbl.Alt(i)*ft2m);
Tbl.rho(i) = rho*0.00194032; % density in slug/ft^3
Tbl.TLapse(i) = TLapse(Tbl.Alt(i)); % thrust lapse
Tbl.KTAS(i) = 565; % true airspeed (kt) (660 mph taken from F18 e/f cruise)
Tbl.KEAS(i) = Tbl.KTAS(i)*sqrt(Tbl.rho(i)/rho_SL); % equivalent airspeed (kt)
Tbl.MACH(i) = Tbl.KTAS(i)/(a*mps2kts); % Mach number

Tbl.CL(i) = (Tbl.WtFrac(i-1)*W_S)/(0.5*Tbl.rho(i)*(Tbl.KTAS(i)*kts2fps)^2); % lift coefficient
Tbl.CD0(i) = 0.008; % drag polar
Tbl.K1(i) = 0; % drag polar
Tbl.K2(i) = 0.08; % drag polar
Tbl.CDR(i) = 0.01; % drag polar
Tbl.CD(i) = Tbl.CD0(i) + Tbl.K1(i)*Tbl.CL(i) + Tbl.K2(i)*Tbl.CL(i)^2 + Tbl.CDR(i); % drag coefficient
Tbl.L_D(i) = Tbl.CL(i)/Tbl.CD(i); % lift-to-drag ratio
Tbl.Drag(i) = Tbl.CD(i) * 0.5* Tbl.rho(i) * (Tbl.KTAS(i)*kts2fps)^2 * (W0/W_S); % drag (lbf)

% solve for throttle setting 
Tbl.Thrust(i) = Tbl.Drag(i); % thrust (lbf)
Tbl.THROT(i) = Tbl.Thrust(i)/(Thrust*Tbl.TLapse(i)); % throttle setting
Tbl.FF(i) = SFC_cruise*Tbl.Thrust(i); % fuel flow (lb/h)
Tbl.dhdt(i) = 0 ; % rate of climb (ft/min)
Tbl.Ps(i) = 0; % specific excess pwr (ft/min)

Tbl.dTime(i) = Tbl.dDist(i)/(Tbl.KTAS(i)/60);% delta time (min)
Tbl.Time(i) = Tbl.dTime(i)+Tbl.Time(i-1); % time (min)
Tbl.FPA(i) = 0; % flightpath angle (deg)
Tbl.GS(i) = Tbl.KTAS(i)*sind(Tbl.FPA(i)); % ground speed (kt)
Tbl.dVdt(i) = 0; % acceleration

Tbl.dFuel(istart9) = 0;
Tbl.dFuel(i) = (Tbl.FF(i)*Tbl.dTime(i))/60; % delta fuel (lb)
Tbl.FuelBurn(i) = sum(Tbl.dFuel(1:i)); % fuel burned (lb)
Tbl.FuelRem(i) = Tbl.FuelRem(i-1) - Tbl.dFuel(i); % fuel remaining (lb)
Tbl.Weight(i) = Tbl.Weight(i-1) - Tbl.dFuel(i); % weight (lb)
Tbl.WtDrop(i) = 0; % dropped weight (lb)
Tbl.WtFrac(i) = Tbl.Weight(i)/Tbl.Weight(i-1);% weight fraction
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 10. descent 2
% 3,000 ft/min

% finish at sea level\
SFC_descent = 0.60;
istart10 = sum(npts(1:9))+1;
iend10 = sum(npts(1:10));
for  i = istart10:iend10

Tbl.Alt(istart10:iend10) = linspace(Tbl.Alt(istart10 - 1),0,npts(10,:)); % altitude (ft)

[T, a, P, rho] = atmosisa(Tbl.Alt(i)*ft2m);
Tbl.rho(i) = rho*0.00194032; % density in slug/ft^3
Tbl.TLapse(i) = TLapse(Tbl.Alt(i)); % thrust lapse
Tbl.MACH(i) = 0.8; % Mach number
Tbl.KTAS(i) = Tbl.MACH(i)*a*mps2kts; % true airspeed (kt)
Tbl.KEAS(i) = Tbl.KTAS(i)*sqrt(Tbl.rho(i)/rho_SL); % equivalent airspeed (kt)

Tbl.THROT(i) = 1; % throttle setting
Tbl.Thrust(i) = Thrust*Tbl.TLapse(i)*Tbl.THROT(i); % thrust (lbf)
Tbl.FF(i) = SFC_descent*Tbl.Thrust(i); % fuel flow (lb/h)
Tbl.dhdt(i) = -3000; % rate of climb (ft/min)
Tbl.Ps(i) = Tbl.dhdt(i); % specific excess pwr (ft/min)

Tbl.CL(i) = (Tbl.WtFrac(i-1)*W_S)/(0.5*Tbl.rho(i)*(Tbl.KTAS(i)*kts2fps)^2); % lift coefficient
Tbl.CD0(i) = 0.008; % drag polar
Tbl.K1(i) = 0; % drag polar
Tbl.K2(i) = 0.08; % drag polar
Tbl.CDR(i) = 0; % drag polar
Tbl.CD(i) = Tbl.CD0(i) + Tbl.K1(i)*Tbl.CL(i) + Tbl.K2(i)*Tbl.CL(i)^2 + Tbl.CDR(i); % drag coefficient
Tbl.L_D(i) = Tbl.CL(i)/Tbl.CD(i); % lift-to-drag ratio
Tbl.Drag(i) = Tbl.CD(i) * 0.5* Tbl.rho(i) * (Tbl.KTAS(i)*kts2fps)^2 * (W0/W_S); % drag (lbf)

Tbl.dTime(i) = (Tbl.Alt(i) - Tbl.Alt(i-1))/Tbl.dhdt(i);% delta time (min)
Tbl.Time(i) = Tbl.dTime(i)+Tbl.Time(i-1); % time (min)
Tbl.GS(i) = Tbl.KTAS(i)*sind(Tbl.FPA(i)); % ground speed (kt)
Tbl.FPA(i) = asind(Tbl.dhdt(i)/(Tbl.KTAS(i)*kts2fps*60)); % flightpath angle (deg)
Tbl.Dist(istart10) = Tbl.Dist(istart10 - 1);
Tbl.Dist(i) = Tbl.GS(i)*Tbl.dTime(i)/60 + Tbl.Dist(istart10 - 1); % distance (NM)
Tbl.dDist(i) =Tbl.Dist(i) - Tbl.Dist(i-1); % delta distance (NM)
Tbl.dVdt(i) = (Tbl.KTAS(i)*kts2fps - (Tbl.KTAS(i-1)*kts2fps))/(Tbl.dTime(i)*60); % acceleration
Tbl.dVdt(istart10) = 0;

Tbl.dFuel(i) = (Tbl.FF(i)*Tbl.dTime(i))/60; % delta fuel (lb)
Tbl.FuelBurn(i) = sum(Tbl.dFuel(1:i)); % fuel burned (lb)
Tbl.FuelRem(i) = Tbl.FuelRem(i-1) - Tbl.dFuel(i); % fuel remaining (lb)
Tbl.Weight(i) = Tbl.Weight(i-1) - Tbl.dFuel(i); % weight (lb)
Tbl.WtDrop(i) = 0; % dropped weight (lb)
Tbl.WtFrac(i) = Tbl.Weight(i)/Tbl.Weight(i-1);% weight fraction

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 11. loiter at sea level
% 20 min duration
% 278 KTAS
% no distance credit
SFC_loiter = 0.75;
istart11 = sum(npts(1:10))+1;
iend11 = sum(npts(1:11));

for  i = istart11:iend11
Tbl.Time(istart11:iend11) = linspace((Tbl.Time(istart11-1)),(Tbl.Time(istart11-1)+20),npts(11,:)); % distance (NM)
Tbl.dTime(istart11) = 0; 
Tbl.dTime(i) = Tbl.Time(i) - Tbl.Time(i-1); % delta distance (NM)
Tbl.Dist(i) = Tbl.Dist(istart);
Tbl.dDist(i) = 0;

Tbl.Alt(i) = 0; % Altitude (ft)
[T, a, P, rho] = atmosisa(Tbl.Alt(i)*ft2m);
Tbl.rho(i) = rho*0.00194032; % density in slug/ft^3
Tbl.TLapse(i) = TLapse(Tbl.Alt(i)); % thrust lapse
Tbl.KTAS(i) = 278; % true airspeed (kt)
Tbl.KEAS(i) = Tbl.KTAS(i)*sqrt(Tbl.rho(i)/rho_SL); % equivalent airspeed (kt)
Tbl.MACH(i) = Tbl.KTAS(i)/(a*mps2kts); % Mach number

Tbl.CL(i) = (Tbl.WtFrac(i-1)*W_S)/(0.5*Tbl.rho(i)*(Tbl.KTAS(i)*kts2fps)^2); % lift coefficient
Tbl.CD0(i) = 0.008; % drag polar
Tbl.K1(i) = 0; % drag polar
Tbl.K2(i) = 0.08; % drag polar
Tbl.CDR(i) = 0; % drag polar
Tbl.CD(i) = Tbl.CD0(i) + Tbl.K1(i)*Tbl.CL(i) + Tbl.K2(i)*Tbl.CL(i)^2 + Tbl.CDR(i); % drag coefficient
Tbl.L_D(i) = Tbl.CL(i)/Tbl.CD(i); % lift-to-drag ratio
Tbl.Drag(i) = Tbl.CD(i) * 0.5* Tbl.rho(i) * (Tbl.KTAS(i)*kts2fps)^2 * (W0/W_S); % drag (lbf)

% solve for throttle setting 
Tbl.Thrust(i) = Tbl.Drag(i); % thrust (lbf)
Tbl.THROT(i) = Tbl.Thrust(i)/(Thrust*Tbl.TLapse(i)); % throttle setting
Tbl.FF(i) = SFC_loiter*Tbl.Thrust(i); % fuel flow (lb/h)
Tbl.dhdt(i) = 0 ; % rate of climb (ft/min)
Tbl.Ps(i) = 0; % specific excess pwr (ft/min)

Tbl.FPA(i) = 0; % flightpath angle (deg)
Tbl.GS(i) = Tbl.KTAS(i)*sind(Tbl.FPA(i)); % ground speed (kt)
Tbl.dVdt(i) = 0; % acceleration

Tbl.dFuel(istart) = 0;
Tbl.dFuel(i) = (Tbl.FF(i)*Tbl.dTime(i))/60; % delta fuel (lb)
Tbl.FuelBurn(i) = sum(Tbl.dFuel(1:i)); % fuel burned (lb)
Tbl.FuelRem(i) = Tbl.FuelRem(i-1) - Tbl.dFuel(i); % fuel remaining (lb)
Tbl.Weight(i) = Tbl.Weight(i-1) - Tbl.dFuel(i); % weight (lb)
Tbl.WtDrop(i) = 0; % dropped weight (lb)
Tbl.WtFrac(i) = Tbl.Weight(i)/Tbl.Weight(i-1);% weight fraction
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 12. landing, taxi-in, shutdown
% 3 min duration

% (T/W)_idle = 0.05;
% no distance credit
SFC_shutdown = 0.6;
istart12 = sum(npts(1:11))+1;
iend12 = sum(npts(1:12));

for  i = istart12:iend12

Tbl.Time(istart12:iend12) = linspace((Tbl.Time(istart12-1)),(Tbl.Time(istart12-1)+ 3),npts(12,:)); % time (min)
Tbl.dTime(istart12) = 0;
Tbl.dTime(i) = Tbl.Time(i) - Tbl.Time(i-1);% delta time (min)
Tbl.Alt(i) = 0; % altitude (ft)
[T, a, P, rho] = atmosisa(Tbl.Alt(i)*ft2m);
Tbl.rho(i) = rho*0.00194032; % density in slug/ft^3
Tbl.KEAS(i) = 0; % equivalent airspeed (kt)
Tbl.KTAS(i) = 0; % true airspeed (kt)
Tbl.MACH(i) = 0; % Mach number
Tbl.GS(i) = 0; % ground speed (kt)
Tbl.FPA(i) = 0; % flightpath angle (deg)
Tbl.Dist(i) = 0; % distance (NM)
Tbl.dDist(i) = 0; % delta distance (NM)
Tbl.dhdt(i) = 0; % rate of climb (ft/min)
Tbl.dVdt(i) = 0; % acceleration
Tbl.CL(i) = 0; % lift coefficient
Tbl.CD0(i) = 0.008; % drag polar
Tbl.K1(i) = 0; % drag polar
Tbl.K2(i) = 0.08; % drag polar
Tbl.CDR(i) = 0; % drag polar
Tbl.CD(i) = 0; % drag coefficient
Tbl.L_D(i) = 0; % lift-to-drag ratio
Tbl.Drag(i) = 0; % drag (lbf)
Tbl.TLapse(i) = 1; % thrust lapse
Tbl.Ps(i) = 0; % specific excess pwr (ft/min)
Tbl.THROT(i) = 1; % throttle setting
Tbl.Thrust(i) = T_W_idle*Tbl.TLapse(i)*W0*Tbl.THROT(i); % thrust (lbf)
Tbl.FF(i) = SFC_shutdown*Tbl.Thrust(i); % fuel flow (lb/h)
Tbl.dFuel(istart12) = 0;
Tbl.dFuel(i) = (Tbl.FF(i)*Tbl.dTime(i))/60; % delta fuel (lb)
Tbl.FuelBurn(i) = sum(Tbl.dFuel(1:i)); % fuel burned (lb)
Tbl.FuelRem(istart) = Tbl.FuelRem(istart - 1);
Tbl.FuelRem(i) = Tbl.FuelRem(i-1) - Tbl.dFuel(i); % fuel remaining (lb)
Tbl.Weight(i) = Tbl.Weight(i-1) - Tbl.dFuel(i); % weight (lb)
Tbl.WtDrop(i) = 0; % dropped weight (lb)
Tbl.WtFrac(i) = Tbl.Weight(i)/Tbl.Weight(i-1);% weight fraction
end





% Any missing info: grab from Exam 2 problem setup
% getting atmosphere properties
% [T_K, a_ms, P_Pa, rho_kgm3] = atmosisa(Alt_m)
% thrust model
% THRUST = [SLS THRUST] * [LAPSE] * [THROT]
% THROT: [0, 1.25]
% 0-1.00: no afterburner
% 1.00 - 1.25: afterburner engaged
% thrust lapse model


% equations

Tbl.EnHt = Tbl.Alt + (Tbl.KTAS*NM2ft).^2/(2*32.17);
Tbl.GS = Tbl.KTAS.*cosd(Tbl.FPA);
disp(Tbl)
A = 2.34;
C = -0.13;
EWF = A*W0^C;
W_crew = 300;
W_payload = 4380; % weight of weapons 
OEW = EWF*W0;
FuelAllow = 0.06*(Tbl.FuelBurn(iend12)); % 6% fuel allowance 
FuelReq = FuelAllow + Tbl.FuelBurn(iend12); 
FuelAvail = W0 - OEW - W_crew - W_payload;
FuelExcess = FuelAvail - FuelReq;
tol = abs(FuelExcess);


W0_calc = FuelReq + OEW + W_crew + W_payload;

W0 = W0_calc;    

end

figure
plot(Tbl.Time, Tbl.FuelBurn)
xlabel("Time (min)")
ylabel("Total Fuel Burn (lb)")

fprintf('Converged Gross Weight is %5.0f lbs .\n', W0)
fprintf('Fuel Required is %5.0f lbs', FuelReq)
