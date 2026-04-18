clearvars;
clc;
close all;

displayResults = false;

valuesPerInput = 2;

WL = linspace(80, 105, valuesPerInput);
AR = linspace(2.5, 5, valuesPerInput);
SA = linspace(25, 55, valuesPerInput);

nCases = valuesPerInput^3;

resultsFlat = cell(nCases,1);
MTOWsFlat = NaN(nCases,1);

parfor idx = 1:nCases
    [i,j,l] = ind2sub([valuesPerInput, valuesPerInput, valuesPerInput], idx);

    aircraft = struct;
    aircraft.constants.wingLoading = WL(i);
    aircraft.wing.AR = AR(j);
    aircraft.wing.taper_ratio = 0.4;
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

results = reshape(resultsFlat, valuesPerInput, valuesPerInput, valuesPerInput);
MTOWs   = reshape(MTOWsFlat,   valuesPerInput, valuesPerInput, valuesPerInput);

%% Plot MTOW as color over (AR, SA, WL)

figure;
hold on;
grid on;
box on;
view(3);

xVals = [];   % AR
yVals = [];   % SA
zVals = [];   % WL
cVals = [];   % MTOW

for i = 1:length(WL)
    for j = 1:length(AR)
        for l = 1:length(SA)

            val = MTOWs(i,j,l);

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
            zVals(end+1) = WL(i);
            cVals(end+1) = val;
        end
    end
end

if ~isempty(cVals)
    scatter3(xVals, yVals, zVals, 100, cVals, "filled");
    cb = colorbar;
    cb.Label.String = "MTOW (lbf)";
else
    text(mean(xlim), mean(ylim), mean(zlim), "No valid data", ...
        "HorizontalAlignment","center");
end

xlabel("Wing Aspect Ratio");
ylabel("Wing LE Sweep Angle [deg]");
zlabel("Wing Loading [lbf/sq.ft]");
title("");

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

%% ============================================================
%% Contour plots: MTOW over (AR, SA), one panel per WL >= 90
% ============================================================

validWLIdx = find(WL >= 90);

nPlots = numel(validWLIdx);
nCols = ceil(sqrt(nPlots));
nRows = ceil(nPlots / nCols);

figure;
tiledlayout(nRows, nCols, "TileSpacing","compact", "Padding","compact");

for p = 1:nPlots
    i = validWLIdx(p);

    nexttile;
    hold on; grid on; box on;

    % Extract MTOW slice for this WL: size = (nAR x nSA)
    M = squeeze(MTOWs(i,:,:));   % AR x SA

    % Ensure numeric and mask invalids
    M = double(M);
    M(~isfinite(M) | ~isreal(M)) = NaN;

    % Build grids (AR along x, SA along y)
    [ARgrid, SAgrid] = meshgrid(AR, SA);

    % Note: M is (AR x SA), but meshgrid is (SA x AR), so transpose M
    Mplot = M';  % now (SA x AR)

    % Plot filled contours
    if any(~isnan(Mplot(:)))
        contourf(ARgrid, SAgrid, Mplot, 12, "LineColor","none");
        colormap(parula);
        cb = colorbar;
        cb.Label.String = "MTOW";
    else
        text(mean(AR), mean(SA), "No valid data", ...
            "HorizontalAlignment","center");
    end

    xlabel("Aspect Ratio, AR");
    ylabel("Sweep Angle, SA [deg]");
    title(sprintf("Wing Loading = %.3g", WL(i)));
end

% Optional: unify color limits across all panels for fair comparison
allVals = MTOWs(WL >= 90, :, :);
allVals = allVals(isfinite(allVals) & isreal(allVals));
if ~isempty(allVals)
    clim([min(allVals), max(allVals)]);
end

%% ============================================================
%% Iso-MTOW contour plots (75k, 85k, 95k) over (AR, SA)
% One panel per WL >= 90
% ============================================================

validWLIdx = find(WL >= 90);

nPlots = numel(validWLIdx);
nCols = ceil(sqrt(nPlots));
nRows = ceil(nPlots / nCols);

isoLevels = [75000:1000:85000];

figure;
tiledlayout(nRows, nCols, "TileSpacing","compact", "Padding","compact");

for p = 1:nPlots
    i = validWLIdx(p);

    nexttile;
    hold on; grid on; box on;

    % Extract MTOW slice (AR x SA)
    M = squeeze(MTOWs(i,:,:));

    % Clean invalid values
    M = double(M);
    M(~isfinite(M) | ~isreal(M)) = NaN;

    % Build grids
    [ARgrid, SAgrid] = meshgrid(AR, SA);
    Mplot = M';  % match meshgrid orientation

    % Only plot if valid data exists
    if any(~isnan(Mplot(:)))
        [C, h] = contour(ARgrid, SAgrid, Mplot, isoLevels, ...
            "LineWidth", 1.8);

        clabel(C, h, "FontSize", 10, "Color", "k");

        legendStrings = arrayfun(@(x) sprintf("%.0f lb", x), isoLevels, ...
            "UniformOutput", false);
        legend(legendStrings, "Location", "best");
    else
        text(mean(AR), mean(SA), "No valid data", ...
            "HorizontalAlignment","center");
    end

    xlabel("Wing Aspect Ratio");
    ylabel("Sweep Angle [deg]");
    title(sprintf("Wing Loading = %.3g lb/sq.ft", WL(i)));
end

%% ============================================================
%% Iso-MTOW contour plot over (AR, SA)
% Only WL = 105
% ============================================================

isoLevels = 75000:1000:85000;

% Find index corresponding to WL = 105 (with tolerance for floating point)
[~, i] = min(abs(WL - 105));

figure;
hold on; grid on; box on;

% Extract MTOW slice (AR x SA)
M = squeeze(MTOWs(i,:,:));

% Clean invalid values
M = double(M);
M(~isfinite(M) | ~isreal(M)) = NaN;

% Build grids
[ARgrid, SAgrid] = meshgrid(AR, SA);
Mplot = M';  % match meshgrid orientation

% Plot contours if valid data exists
if any(~isnan(Mplot(:)))
    [C, h] = contour(ARgrid, SAgrid, Mplot, isoLevels, ...
        "LineWidth", 1.6);

    % Label contours with units
    clabel(C, h, ...
        "FontSize", 10, ...
        "Color", "k")%, ...
        %"LabelFormat", "%.0f lbf");
else
    text(mean(AR), mean(SA), "No valid data", ...
        "HorizontalAlignment","center");
end

xlabel("Aspect Ratio, AR");
ylabel("Sweep Angle, SA [deg]");
title(sprintf("Iso-MTOW Contours at WL = %.0f lbf/ft^2", WL(i)));

%% ============================================================
%% Iso-MTOW contour plot over (AR, SA)
% Only WL = 105
% ============================================================

isoLevels = 72000:1000:88000;

% Find index corresponding to WL = 105
[~, i] = min(abs(WL - 105));

figure;
hold on;
grid on;
box on;

% Extract MTOW slice (AR x SA)
M = squeeze(MTOWs(i,:,:));

% Clean invalid values
M = double(M);
M(~isfinite(M) | ~isreal(M)) = NaN;

% Build grids
[ARgrid, SAgrid] = meshgrid(AR, SA);
Mplot = M';   % match meshgrid orientation

if any(~isnan(Mplot(:)))
    [C, h] = contour(ARgrid, SAgrid, Mplot, isoLevels, "LineWidth", 2.8);

    % Basic contour labels
    clabel(C, h, "FontSize", 12, "Color", "k", "LabelSpacing", 600);

    % Legend with units
    legend(h(1), "MTOW contours [lbf]", "Location", "bestoutside");
else
    text(mean(AR), mean(SA), "No valid data", ...
        "HorizontalAlignment", "center");
end

xlabel("Wing Aspect Ratio");
ylabel("Wing LE Sweep Angle [deg]");
% title(sprintf("Iso-MTOW Contours at WL = %.0f lbf/ft^2", WL(i)));