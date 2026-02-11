%% runMission_master.m
clc; clear;

cfg = getConfig(); %#ok<NASGU>

SegNames = {'SWT','TKO','CLIMB 1','CR OBD','DESC 1','COMBAT','WP FIRE','CLIMB 2', ...
            'CR IBD','DESC 2','LTR2','LTS'};

npts = [ ...
    3;  % 1 SWT
    3;  % 2 TKO
    10; % 3 CLIMB 1
    10; % 4 CR OBD
    10; % 5 DESC 1
    5;  % 6 COMBAT
    1;  % 7 WP FIRE
    10;     10; % 8 CLIMB 2
    10; % 9 CR IBD
    10; % 10 DESC 2
    10; % 11 LTR2
    3;  % 12 LTS
];

W0 = 63738;
FuelReq = 16440;
tol = 2;

while tol > 1
    [Tbl, FuelReq_new, W0_next] = runMissionOnce(W0, FuelReq, SegNames, npts);

    tol = abs(W0_next - W0);

    W0 = W0_next;
    FuelReq = FuelReq_new;
end

figure
plot(Tbl.Time, Tbl.FuelBurn)
xlabel("Time (min)")
ylabel("Total Fuel Burn (lb)")

disp(Tbl)
fprintf('Converged Gross Weight is %5.0f lbs .\n', W0)
fprintf('Fuel Required is %5.0f lbs\n', FuelReq)
