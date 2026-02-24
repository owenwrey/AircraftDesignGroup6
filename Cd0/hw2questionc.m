%% Problem 3(d): Euler integration for xdot = -10x, x(0)=1 on [0,2]

clear; clc; close all;

% given
t0 = 0;            % initial time [s]
tf = 2;            % final time [s]
x0 = 1;            % initial condition
dts = [0.25 0.2 0.15 0.02];   % time steps to test

% ODE
f = @(t,x) -10*x;

% exact solution for comparison (smooth time grid)
t_exact = linspace(t0, tf, 2001);
x_exact = exp(-10*t_exact);

% loop over each dt and plot in a separate figure
for i = 1:numel(dts)
    dt = dts(i);

    % time grid for this dt
    t = t0:dt:tf;
    N = numel(t);

    % preallocate solution
    x = zeros(1, N);
    x(1) = x0;

    % Forward Euler integration
    for k = 1:N-1
        x(k+1) = x(k) + dt * f(t(k), x(k));
    end

    % plot this dt in its own figure
    figure; grid on; hold on;
    plot(t_exact, x_exact, 'LineWidth', 2);              % exact
    plot(t, x, 'o-', 'LineWidth', 1.2, 'MarkerSize', 5); % euler

    xlabel('t [s]');
    ylabel('x(t)');
    title(sprintf('Forward Euler vs Exact: \\Deltat = %.2f', dt));
    legend('Exact: e^{-10t}', sprintf('Euler: \\Deltat = %.2f', dt), 'Location', 'best');
end
