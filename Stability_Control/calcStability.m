clearvars;
close all;
% clc;

tbl = readtable("dataXLSX1.xlsx");

tbl.CD = tbl.CDi + tbl.CD0;

tbl.Cl = tbl.CMx;
tbl.Cm = tbl.CMy;
tbl.Cn = tbl.CMz;

s = struct;

% make constant beta sub structs
betaNames = {'b0','b2','b4','b6'};
betas = [0,2,4,6];

for i = 1:numel(betas)
    s.(betaNames{i}).tbl = tbl(tbl.beta == betas(i), :);
end

% make constant AOA sub structs
alphaNames = {'an10','an5','a0','a5','a10','a15','a20','a25'};
alphas = -10:5:25;

for i = 1:numel(alphas)
    s.(alphaNames{i}).tbl = tbl(tbl.AOA == alphas(i), :);
end

%% CLalpha
figure;
plot(s.b0.tbl.AOA, s.b0.tbl.CL,'LineWidth',1.75); grid on;
xlabel('Angle of attack, $\alpha$ (deg)','Interpreter','latex')
ylabel('Lift Coefficient, $C_L$','Interpreter','latex')

%% CLbeta

%% CDalpha

%% CMalpha
figure;
plot(s.b2.tbl.AOA, s.b2.tbl.Cm,'LineWidth',1.75); grid on;
xlabel('Angle of attack, $\alpha$ (deg)','Interpreter','latex')
ylabel('Pitching Moment Coefficient, $C_m$','Interpreter','latex')


%% CYbeta
figure;
plot(s.a5.tbl.beta, s.a5.tbl.Cy,'LineWidth',1.75,'Color',[0, .55, 0]); grid on;
xlabel('Sideslip Angle, $\beta$ (deg)','Interpreter','latex')
ylabel('Side Force Coefficient, $C_y$','Interpreter','latex')


%% Clbeta
figure;
plot(s.a5.tbl.beta, s.a5.tbl.Cl,'LineWidth',1.75,'Color',[0, .55, 0]); grid on;
xlabel('Sideslip Angle, $\beta$ (deg)','Interpreter','latex')
ylabel('Rolling Moment Coefficient, $C_l$','Interpreter','latex')


%% Cnbeta
figure;
plot(s.a5.tbl.beta, s.a5.tbl.Cn,'LineWidth',1.75,'Color',[0, .55, 0]); grid on;
xlabel('Sideslip Angle, $\beta$ (deg)','Interpreter','latex')
ylabel('Yaw Moment Coefficient, $C_n$','Interpreter','latex')


%% Alpha plots
% CL vs alpha
figure;
tiledlayout(1,2,"TileSpacing","compact","Padding","compact");
nexttile;
plot(s.b0.tbl.AOA, s.b0.tbl.CL,'LineWidth',1.75);
grid on;
xlabel('Angle of attack, $\alpha$ (deg)','Interpreter','latex')
ylabel('$C_L$','Interpreter','latex')
title('$C_L$ vs $\alpha$','Interpreter','latex')


% Cm vs alpha
nexttile;
plot(s.b2.tbl.AOA, s.b2.tbl.Cm,'LineWidth',1.75);
grid on;
xlabel('Angle of attack, $\alpha$ (deg)','Interpreter','latex')
ylabel('$C_m$','Interpreter','latex')
title('$C_m$ vs $\alpha$','Interpreter','latex')

%% Beta plots
figure;
tiledlayout(1,3,"TileSpacing","compact","Padding","compact");

% Cy vs beta
nexttile;
plot(s.a5.tbl.beta, s.a5.tbl.Cy,'LineWidth',1.75,'Color',[0, .55, 0]);
grid on;
xlabel('Sideslip Angle, $\beta$ (deg)','Interpreter','latex')
ylabel('$C_y$','Interpreter','latex')
title('$C_y$ vs $\beta$','Interpreter','latex')

% Cl vs beta
nexttile;
plot(s.a5.tbl.beta, s.a5.tbl.Cl,'LineWidth',1.75,'Color',[0, .55, 0]);
grid on;
xlabel('Sideslip Angle, $\beta$ (deg)','Interpreter','latex')
ylabel('$C_l$','Interpreter','latex')
title('$C_l$ vs $\beta$','Interpreter','latex')

% Cn vs beta
nexttile;
plot(s.a5.tbl.beta, s.a5.tbl.Cn,'LineWidth',1.75,'Color',[0, .55, 0]);
grid on;
xlabel('Sideslip Angle, $\beta$ (deg)','Interpreter','latex')
ylabel('$C_n$','Interpreter','latex')
title('$C_n$ vs $\beta$','Interpreter','latex')

%% PRint derivs
CLalpha = ( (s.b0.tbl.CL(end/2))-(s.b0.tbl.CL(end/2-1)) ) / ((s.b0.tbl.AOA(end/2))-(s.b0.tbl.AOA(end/2-1))) *180/pi;
fprintf('CLalpha = %.3f /rad\n', CLalpha);

CMalpha = ( (s.b2.tbl.Cm(end/2+1))-(s.b2.tbl.Cm(end/2)) ) / ((s.b0.tbl.AOA(end/2))-(s.b0.tbl.AOA(end/2-1))) *180/pi;
fprintf('CMalpha = %.3f /rad\n', CMalpha);

CYbeta = ( (s.a5.tbl.Cy(end/2))-(s.a5.tbl.Cy(end/2-1)) ) / ((s.a5.tbl.beta(3))-(s.a5.tbl.beta(2))) *180/pi;
fprintf('CYbeta = %.3f /rad\n', CYbeta);

Clbeta = ( (s.a0.tbl.Cl(end/2))-(s.a0.tbl.Cl(end/2-1)) ) / ((s.a0.tbl.beta(3))-(s.a0.tbl.beta(2))) *180/pi;
fprintf('Clbeta = %.3f /rad\n', Clbeta);

Cnbeta = ( (s.a0.tbl.Cn(end/2))-(s.a0.tbl.Cn(end/2-1)) ) / ((s.a0.tbl.beta(3))-(s.a0.tbl.beta(2))) *180/pi;
fprintf('Cnbeta = %.3f /rad\n', Cnbeta);