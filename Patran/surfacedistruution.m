clear; clc;

%% Load data
data = readmatrix("FlightStreamPressure.csv");

x = data(:,1);
y = data(:,2);
z = data(:,3);
Cp = data(:,4);

%% === 1. Pick span location (example: mid-span) ===
y_target = mean(y);

tol = 0.01 * (max(y) - min(y)); % tolerance band

idx = abs(y - y_target) < tol;

x_slice = x(idx);
z_slice = z(idx);
Cp_slice = Cp(idx);

%% === 2. Normalize chord ===
x_min = min(x_slice);
x_max = max(x_slice);

x_norm = (x_slice - x_min) / (x_max - x_min);

%% === 3. Separate upper surface ===
% assumes z is vertical direction
upper_idx = z_slice >= mean(z_slice);

x_upper = x_norm(upper_idx);
Cp_upper = Cp_slice(upper_idx);

%% === 4. Sort for clean curve ===
[x_upper, sort_idx] = sort(x_upper);
Cp_upper = Cp_upper(sort_idx);

%% === 5. Export ===
UpperSurfacePressureDist = [x_upper, Cp_upper];

writematrix(UpperSurfacePressureDist, "UpperSurfacePressureDist.csv");

disp("Upper surface pressure distribution exported.");