clear all; clc;

aircraft = struct;

% Ensure required sub-structs exist
if ~isfield(aircraft, 'constants') || ~isstruct(aircraft.constants)
    aircraft.constants = struct();
end

if ~isfield(aircraft, 'weight') || ~isstruct(aircraft.weight)
    aircraft.weight = struct();
end

% Only set defaults if they do not already exist
if ~isfield(aircraft.constants, 'wingLoading')
    aircraft.constants.wingLoading = 115;
    warning("Used hardcoded wingLoading")
end

if ~isfield(aircraft.weight, 'tolerance')
    aircraft.weight.tolerance = 0.06;
    warning("Used hardcoded weight.tolerance")
end


% aircraft = A2A(aircraft);
aircraft = Strike(aircraft);