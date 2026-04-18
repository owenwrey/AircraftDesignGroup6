clear; clc;

%% === 1. Load data ===
data = readmatrix("FA37_pressureProfile_v1.csv");

% Correct column mapping (based on your data)
y_span = data(:,1);   % span
x_chord = data(:,2);  % chord
z = data(:,3);        % thickness (up/down)
Cp = data(:,4);       % pressure

%% === 2. Pick span slice (mid-span) ===
y_target = mean(y_span);

tol = 0.05 * (max(y_span) - min(y_span));  % increase tolerance for safety

idx = abs(y_span - y_target) < tol;

x_slice = x_chord(idx);
z_slice = z(idx);
Cp_slice = Cp(idx);

%% === 3. Normalize chord (0 → 1) ===
x_min = min(x_slice);
x_max = max(x_slice);

x_norm = (x_slice - x_min) / (x_max - x_min);

%% === 4. Split upper and lower surfaces ===
z_mid = (max(z_slice) + min(z_slice)) / 2;

upper_idx = z_slice > z_mid;
lower_idx = z_slice <= z_mid;

x_upper = x_norm(upper_idx);
Cp_upper = Cp_slice(upper_idx);

x_lower = x_norm(lower_idx);
Cp_lower = Cp_slice(lower_idx);

%% === 5. Sort (important for clean curve) ===
[x_upper, i1] = sort(x_upper);
Cp_upper = Cp_upper(i1);

[x_lower, i2] = sort(x_lower);
Cp_lower = Cp_lower(i2);

%% === 6. Export CSV files ===
writematrix([x_upper, Cp_upper], "UpperSurfacePressureDist.csv");
writematrix([x_lower, Cp_lower], "LowerSurfacePressureDist.csv");

disp("Upper and Lower surface CSV files created.");

%% === 7. Plot to verify ===
figure;
scatter(x_upper, Cp_upper, 'r', 'filled'); hold on;
scatter(x_lower, Cp_lower, 'b', 'filled');
xlabel('x/c'); ylabel('Cp');
legend('Upper Surface','Lower Surface');
title('Pressure Distribution');
grid on;