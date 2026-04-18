clear; clc;

%% Load data
data = readmatrix("FlightStreamPressure.csv");

x = data(:,1);
y = data(:,2);
z = data(:,3);
Cp = data(:,4);

%% === 1. Pick span slice ===
y_target = mean(y);
tol = 0.01 * (max(y) - min(y));

idx = abs(y - y_target) < tol;

x_slice = x(idx);
z_slice = z(idx);
Cp_slice = Cp(idx);

%% === 2. Normalize chord ===
x_min = min(x_slice);
x_max = max(x_slice);
x_norm = (x_slice - x_min) / (x_max - x_min);

%% === 3. Split upper/lower using z ===
z_mid = (max(z_slice) + min(z_slice)) / 2;

upper_idx = z_slice > z_mid;
lower_idx = z_slice <= z_mid;

x_upper = x_norm(upper_idx);
Cp_upper = Cp_slice(upper_idx);

x_lower = x_norm(lower_idx);
Cp_lower = Cp_slice(lower_idx);

%% === 4. Sort ===
[x_upper, i1] = sort(x_upper);
Cp_upper = Cp_upper(i1);

[x_lower, i2] = sort(x_lower);
Cp_lower = Cp_lower(i2);

%% === 5. Export ===
writematrix([x_upper, Cp_upper], "UpperSurfacePressureDist.csv");
writematrix([x_lower, Cp_lower], "LowerSurfacePressureDist.csv");

disp("Upper and lower distributions created.");