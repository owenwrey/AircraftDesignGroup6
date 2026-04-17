clc;
clear;
close all;

%% ============================================================
%  NSGA-II / Multi-objective GA tutorial in MATLAB
%  Using MATLAB built-in multi-objective GA: gamultiobj
%  ============================================================

%% ---------------- DESIGN VARIABLES ----------------
% x(1) = Wing loading
% x(2) = Aspect ratio
% ....
% ...
% x(n) = .......

% x(1) = WingLoading;
% x(2) = AspectRatio;
% x(3) = TaperRatio;
% x(4) = sweepAngle;

     
% lb refers to lower bound of the design variable and ub refers to the
% upper bound of the design variables. In this example the lower bound of
% wing loading in 10 units and that for aspect ratio is 10. The upper bound
% of wing loading is 30 units and that for aspect ratio is 15. lb and ub
% should be a row vector of 1xnvars
lb = [80, 2.5, 0.2, 25];    
ub = [105, 6, 0.45, 50];

nvars = length(lb); % It should be equal to the number of design variables

% Population size is the total number of designs that the optimizer will
% evaluate in each generation. A rule of thumb is to make populationSize =
% 5xnvars.
populationSize = max(20, 5 * nvars);

% Maintain the crossover fraction to 0.95. Please do not change this
crossoverFraction = 0.95;

%% ---------------- OPTIONS ----------------
options = optimoptions('gamultiobj', ...
    'FunctionTolerance', 1e-3, ...
    'MaxStallGenerations', 20, ...
    'UseParallel', true, ...
    'PopulationSize',    populationSize, ...
    'CrossoverFraction', crossoverFraction, ...
    'MaxGenerations',    100, ...
    'FunctionTolerance', 1e-4, ...
    'MutationFcn',       @mutationadaptfeasible, ...
    'Display',           'iter', ...
    'OutputFcn',         @outputfunction ...
    );  % <-- your custom output function added here

%% ---------------- RUN OPTIMIZATION ----------------
[xPareto, fPareto, exitflag, output, population, scores] = ...
    gamultiobj(@objectiveFunction, nvars, [], [], [], [], lb, ub, [], options);

%% ============================================================
% OBJECTIVE FUNCTION
%% ============================================================
function f = objectiveFunction(x)

    persistent errorLog

    if isempty(errorLog)
        errorLog = {};
    end

    try
        aircraft = SizingFrameworkApr09(x);
    
        % Extract the necessary metrics that you want as objective functions to
        % minimize or maximize
        obj1 = aircraft.weight.total;
        obj2 = aircraft.constants.EWF;
        obj3 = aircraft.aero.cd0_strike_interp(1.2, 30000);
        
        % If the objective function is to minimize (like MTOM), do not change
        % the sign. But, if the objective function is to maximize (like payload
        % mass fraction), change the sign to -obj2. If the second objective function 
        % was to minimize fuel weight, it would just be 0bj2, no change in sign. 
        % You can add one more objective function if there is any. This setup can handle only 3
        % objective functions
        
        f = [obj1, obj2, obj3];
    catch ERRMSG
        % Append error info
        errorLog{end+1} = struct( ...
            'x', x, ...
            'message', ERRMSG.message, ...
            'identifier', ERRMSG.identifier, ...
            'stack', ERRMSG.stack);

        f = [1e12, 1000, 1000];

        % warning("Failure in objectiveFunction: %s", ERRMSG.message)
    end
end


