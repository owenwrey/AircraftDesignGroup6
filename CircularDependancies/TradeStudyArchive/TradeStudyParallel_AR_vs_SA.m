clearvars;
clc;
close all;

displayResults = false;

valuesPerInput = 10;

WL = 102;
AR = linspace(2.5, 5, valuesPerInput);
SA = linspace(10, 45, valuesPerInput);

nCases = valuesPerInput^2;

resultsFlat = cell(nCases,1);
MTOWsFlat = NaN(nCases,1);

parfor idx = 1:nCases
    [j,l] = ind2sub([valuesPerInput, valuesPerInput], idx);

    aircraft = struct;
    aircraft.constants.wingLoading = WL;
    aircraft.wing.AR = AR(j);
    aircraft.wing.taper_ratio = 0.3;
    aircraft.wing.sweep = SA(l);

    try
        aircraft = Main_forTS(aircraft);
        resultsFlat{idx} = aircraft;
        MTOWsFlat(idx) = aircraft.weight.total;
    catch ERRMSG
        resultsFlat{idx} = ERRMSG.message;
        MTOWsFlat(idx) = NaN;
    end
end

results = reshape(resultsFlat, valuesPerInput, valuesPerInput);
MTOWs   = reshape(MTOWsFlat,   valuesPerInput, valuesPerInput);

%% Plot MTOW as color over (AR, SA)

figure;
hold on;
grid on;
box on;

xVals = [];   % AR
yVals = [];   % SA
cVals = [];   % MTOW

for j = 1:length(AR)
    for l = 1:length(SA)

        val = MTOWs(j,l);

        % Skip invalid values
        if isempty(val)
            continue
        end
        if ~isnumeric(val) || ~isscalar(val)
            continue
        end
        if ~isreal(val) || isnan(val)
            continue
        end

        xVals(end+1) = AR(j);
        yVals(end+1) = SA(l);
        cVals(end+1) = val;
    end
end

if ~isempty(cVals)
    scatter(xVals, yVals, 100, cVals, "filled");
    cb = colorbar;
    cb.Label.String = "MTOW (lbf)";
else
    text(mean(xlim), mean(ylim), "No valid data", ...
        "HorizontalAlignment","center");
end

xlabel("Wing Aspect Ratio");
ylabel("Wing LE Sweep Angle [deg]");
title(sprintf("WL = %.0f lbf/sq.ft", WL));

%% Count error types
resultsFlat = results(:);

isText = cellfun(@(x) ischar(x) || isstring(x), resultsFlat);
textEntries = string(resultsFlat(isText));
textEntries = textEntries(textEntries ~= "");

[uniqueMsgs, ~, idx] = unique(textEntries);
counts = accumarray(idx, 1);

[counts, order] = sort(counts, 'descend');
uniqueMsgs = uniqueMsgs(order);

for n = 1:length(uniqueMsgs)
    fprintf("%5d  :  %s\n", counts(n), uniqueMsgs(n));
end

%% Iso-MTOW contour plot over (AR, SA)

isoLevels = 72000:500:88000;

figure;
hold on;
grid on;
box on;
axis square;

% MTOW slice is now just AR x SA
M = MTOWs;

% Clean invalid values
M = double(M);
M(~isfinite(M) | ~isreal(M)) = NaN;

% Build grids
[ARgrid, SAgrid] = meshgrid(AR, SA);
Mplot = M';   % match meshgrid orientation

if any(~isnan(Mplot(:)))
    [C, h] = contour(ARgrid, SAgrid, Mplot, isoLevels, "LineWidth", 2.8);
    clabel(C, h, "FontSize", 12, "Color", "k", "LabelSpacing", 600);
else
    text(mean(AR), mean(SA), "No valid data", ...
        "HorizontalAlignment", "center");
end

xlabel("Wing Aspect Ratio", Interpreter="latex");
ylabel("Wing Quarter-chord Sweep Angle, [deg]", Interpreter="latex");

% Plot a single design point (manual input)

% ---- User inputs ----
WL_pt   = 102;      % lbf/ft^2
AR_pt   = 3.25;
SA_pt   = 25;       % deg
MTOW_pt = 74653;    % lbf

% ---- Plot ----
hPt = plot(AR_pt, SA_pt, 'r*', ...
    'DisplayName', 'Selected Design Point', ...
    'MarkerSize', 14, ...
    'MarkerFaceColor', 'r', ...
    'LineWidth', 1.5);

% Optional: check consistency with current WL
if abs(WL_pt - WL) > 1e-6
    warning('Plotted point WL does not match current analysis WL');
end

% Convert data coordinates -> normalized figure coordinates
ax = gca;
pt = [AR_pt, SA_pt];
pt_norm = ax.Position(1:2) + ...
    [(pt(1)-ax.XLim(1))/diff(ax.XLim)*ax.Position(3), ...
     (pt(2)-ax.YLim(1))/diff(ax.YLim)*ax.Position(4)];

% Create label string
labelStr = { ...
    sprintf('$W/S = %.0f\\ \\mathrm{lb_f/ft^2}$', WL_pt), ...
    sprintf('$AR = %.2f$', AR_pt), ...
    sprintf('$\\Lambda_{c/4} = %.1f^{\\circ}$', SA_pt), ...
    sprintf('$MTOW = %.0f\\ \\mathrm{lb_f}$', MTOW_pt) ...
    };

% Add annotation box
annotation('textbox', ...
    [pt_norm(1)+0.055, pt_norm(2)+0.05, 0.15, 0.1], ...
    'String', labelStr, ...
    'FitBoxToText', 'on', ...
    'BackgroundColor', 'w', ...
    'EdgeColor', 'k', ...
    'FontSize', 12, ...
    'Interpreter', 'latex');

% Explicit legend entries
if any(~isnan(Mplot(:)))
    legend([h(1), hPt], {"MTOW contours [lbf]", "Selected Design Point"}, ...
        "Location", "northwest");
else
    legend(hPt, "Selected Design Point", "Location", "northeast");
end