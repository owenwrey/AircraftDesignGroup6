clearvars;
clc;

displayResults = false;

valuesPerInput = 6;

WL = linspace(80, 105, valuesPerInput);
AR = linspace(2.5, 5, valuesPerInput);
TR = linspace(0.2, 0.55, valuesPerInput);
SA = linspace(25, 55, valuesPerInput);

nCases = valuesPerInput^4;

resultsFlat = cell(nCases,1);
MTOWsFlat = NaN(nCases,1);

parfor idx = 1:nCases
    [i,j,k,l] = ind2sub([valuesPerInput, valuesPerInput, valuesPerInput, valuesPerInput], idx);

    aircraft = struct;
    aircraft.constants.wingLoading = WL(i);
    aircraft.wing.AR = AR(j);
    aircraft.wing.taper_ratio = TR(k);
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

results = reshape(resultsFlat, valuesPerInput, valuesPerInput, valuesPerInput, valuesPerInput);
MTOWs   = reshape(MTOWsFlat,   valuesPerInput, valuesPerInput, valuesPerInput, valuesPerInput);

%% Plot MTOW as color over (AR, TR, SA), one figure panel per WL

figure;
tiledlayout(2,3, "TileSpacing","compact", "Padding","compact");

for i = 1:length(WL)
    nexttile;
    hold on;
    grid on;
    box on;
    view(3);

    xVals = [];   % AR
    yVals = [];   % TR
    zVals = [];   % SA
    cVals = [];   % MTOW

    for j = 1:length(AR)
        for k = 1:length(TR)
            for l = 1:length(SA)

                % Get MTOW value robustly for either cell or numeric array
                if iscell(MTOWs)
                    val = MTOWs{i,j,k,l};
                else
                    val = MTOWs(i,j,k,l);
                end

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
                yVals(end+1) = TR(k);
                zVals(end+1) = SA(l);
                cVals(end+1) = val;
            end
        end
    end

    if ~isempty(cVals)
        scatter3(xVals, yVals, zVals, 100, cVals, "filled");
        cb = colorbar;
        cb.Label.String = "MTOW";
    else
        text(mean(xlim), mean(ylim), mean(zlim), "No valid data", ...
            "HorizontalAlignment","center");
    end

    xlabel("Wing Aspect Ratio");
    ylabel("Wing Taper Ratio");
    zlabel("Wing LE Sweep Angle [deg]");
    title(sprintf("Wing Loading = %.3g lbf/sq.ft", WL(i)));
end

%% Count error types
% Flatten
resultsFlat = results(:);

% Extract only text entries (string or char)
isText = cellfun(@(x) ischar(x) || isstring(x), resultsFlat);

textEntries = string(resultsFlat(isText));  % normalize to string

% Optional: remove empty strings
textEntries = textEntries(textEntries ~= "");

% Count occurrences
[uniqueMsgs, ~, idx] = unique(textEntries);
counts = accumarray(idx, 1);

% Sort by frequency (descending)
[counts, order] = sort(counts, 'descend');
uniqueMsgs = uniqueMsgs(order);

% Display
for n = 1:length(uniqueMsgs)
    fprintf("%5d  :  %s\n", counts(n), uniqueMsgs(n));
end

%% MTOW vs TR at WL = 105, colored by AR and SA

[~, iWL] = min(abs(WL - 105));
M = squeeze(MTOWs(iWL,:,:,:)); % [AR, TR, SA]

figure;
tiledlayout(1,2, "TileSpacing","compact", "Padding","compact");

% Two Plots that show that TR doesn't change much
% -------- Plot 1: colored by AR --------
nexttile;
axis square;
hold on; 
grid on; 
box on;

cmapAR = nebula(length(AR));

for j = 1:length(AR)
    xVals = [];
    yVals = [];

    for k = 1:length(TR)
        for l = 1:length(SA)
            val = M(j,k,l);

            if ~isfinite(val)
                continue
            end

            xVals(end+1) = TR(k);
            yVals(end+1) = val;
        end
    end

    scatter(xVals, yVals, 80, ...
        'MarkerFaceColor', cmapAR(j,:), ...
        'MarkerEdgeColor','none',...
        'DisplayName', sprintf("AR = %.2f", AR(j)));
end

xlim([min(TR) - 0.05*(max(TR)-min(TR)), max(TR) + 0.05*(max(TR)-min(TR))])

xlabel("Wing Taper Ratio", Interpreter="latex");
ylabel("MTOW [lbf]", Interpreter="latex");
title(sprintf("MTOW vs Taper Ratio (Colored by AR, W/S = %.3g lb/ft^2)", WL(iWL)));

legend("Location","northoutside","Orientation","horizontal");


% -------- Plot 2: colored by SA --------
nexttile;
axis square
hold on; 
grid on; 
box on;

cmapSA = parula(length(SA));  % different scheme

for l = 1:length(SA)
    xVals = [];
    yVals = [];

    for j = 1:length(AR)
        for k = 1:length(TR)
            val = M(j,k,l);

            if ~isfinite(val)
                continue
            end

            xVals(end+1) = TR(k);
            yVals(end+1) = val;
        end
    end

    scatter(xVals, yVals, 80, ...
        'MarkerFaceColor', cmapSA(l,:), ...
        'MarkerEdgeColor','none',...
        'DisplayName', sprintf("SA = %.1f°", SA(l)));
end

xlim([min(TR) - 0.05*(max(TR)-min(TR)), max(TR) + 0.05*(max(TR)-min(TR))])

xlabel("Wing Taper Ratio", Interpreter="latex");
ylabel("MTOW [lbf]", Interpreter="latex");
title(sprintf("MTOW vs Taper Ratio (Colored by Sweep Angle, W/S = %.3g lb/ft^2)", WL(iWL)));

legend("Location","northoutside","Orientation","horizontal");

%% MTOW vs WL at fixed TR, colored by AR and SA

[~, iTR] = min(abs(TR - 0.34));
M = squeeze(MTOWs(:,:,iTR,:)); % [WL, AR, SA]

figure;
tiledlayout(1,2, "TileSpacing","compact", "Padding","compact");

% Two Plots that show how WL changes MTOW
% -------- Plot 1: colored by AR --------
nexttile;
axis square;
hold on;
grid on;
box on;

cmapAR = nebula(length(AR));

for j = 1:length(AR)
    xVals = [];
    yVals = [];

    for i = 1:length(WL)
        for l = 1:length(SA)
            val = M(i,j,l);

            if ~isfinite(val)
                continue
            end

            xVals(end+1) = WL(i);
            yVals(end+1) = val;
        end
    end

    scatter(xVals, yVals, 80, ...
        'MarkerFaceColor', cmapAR(j,:), ...
        'MarkerEdgeColor', 'none', ...
        'DisplayName', sprintf("AR = %.2f", AR(j)));
end

xlim([min(WL) - 0.05*(max(WL)-min(WL)), max(WL) + 0.05*(max(WL)-min(WL))])

xlabel("Wing Loading, W/S [lb/ft^2]", Interpreter="latex");
ylabel("MTOW [lbf]", Interpreter="latex");
title(sprintf("MTOW vs Wing Loading (Colored by AR, TR = %.3g)", TR(iTR)));

legend("Location","northoutside","Orientation","horizontal");


% -------- Plot 2: colored by SA --------
nexttile;
axis square;
hold on;
grid on;
box on;

cmapSA = parula(length(SA));

for l = 1:length(SA)
    xVals = [];
    yVals = [];

    for i = 1:length(WL)
        for j = 1:length(AR)
            val = M(i,j,l);

            if ~isfinite(val)
                continue
            end

            xVals(end+1) = WL(i);
            yVals(end+1) = val;
        end
    end

    scatter(xVals, yVals, 80, ...
        'MarkerFaceColor', cmapSA(l,:), ...
        'MarkerEdgeColor', 'none', ...
        'DisplayName', sprintf("SA = %.1f°", SA(l)));
end

xlim([min(WL) - 0.05*(max(WL)-min(WL)), max(WL) + 0.05*(max(WL)-min(WL))])

xlabel("Wing Loading, W/S [lb/ft^2]", Interpreter="latex");
ylabel("MTOW [lbf]", Interpreter="latex");
title(sprintf("MTOW vs Wing Loading (Colored by Sweep Angle, TR = %.3g)", TR(iTR)));

legend("Location","northoutside","Orientation","horizontal");

