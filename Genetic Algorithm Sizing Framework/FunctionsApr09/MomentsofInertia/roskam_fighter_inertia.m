function aircraft = inertia(aircraft)
% roskam_fighter_inertia
% Selects a fighter analog from Roskam Table B9a and computes
% moments of inertia at MTOW and empty weight
%
% Inputs:
%   MTOW_lb   - maximum takeoff weight [lb]
%   Wempty_lb - empty weight [lb]
%   b_ft      - wingspan [ft]
%   L_ft      - total length [ft]
%
% Output:
%   inertia   - struct with chosen reference aircraft, nondimensional radii,
%               and inertia values at MTOW and empty weight
%
% Units:
%   inertia outputs are in slug-ft^2

    g = 32.2; % ft/s^2

    % -----------------------------
    % Roskam Table B9a fighter jet data
    % -----------------------------
    db = struct( ...
        'name', { ...
            'McD F2H-1','McD F3H-2N','McD F-101A','VS Attacker','DH Vampire 20', ...
            'Gloster Meteor II','Lockheed F-80A','Lockheed F-94B','Lockheed F-104G', ...
            'NAA F-86A','NAA FJ-3','NAA F-100D','Vought XF8U-1','Vought F8U-3', ...
            'GD XF-91','GD TF-102A','GD F-106B','Northrop F-89D','Republic RF-84F', ...
            'Republic F-105D','Grumman F9F-8','Grumman XF10F-1','Grumman F11F-1' ...
        }, ...
        'GW_lb', { ...
            14413,26878,36969,10450,10891,11100,11940,13650,20900,13900,16883,29800,21300,30600,18600,32859,36834,38000,19000,34058,16744,26160,16500 ...
        }, ...
        'b_ft', { ...
            41.6,35.3,39.7,36.9,40.0,43.0,38.9,37.5,21.9,37.1,37.0,38.0,35.7,40.0,31.3,38.1,38.3,58.0,33.6,35.0,34.5,36.8,31.6 ...
        }, ...
        'L_ft', { ...
            40.2,58.8,67.4,37.3,30.1,41.4,34.3,40.1,54.8,37.5,37.6,47.0,54.4,58.9,43.3,63.2,70.7,54.0,47.5,64.4,41.9,57.8,40.8 ...
        }, ...
        'Rx_bar', { ...
            0.230,0.252,0.209,0.244,0.286,0.286,0.286,0.284,0.224,0.266,0.281,0.252,0.225,0.225,0.323,0.295,0.247,0.440,0.310,0.231,0.248,0.251,0.221 ...
        }, ...
        'Ry_bar', { ...
            0.359,0.107,0.329,0.328,0.318,0.330,0.356,0.396,0.392,0.346,0.352,0.376,0.404,0.375,0.424,0.386,0.379,0.304,0.308,0.425,0.374,0.323,0.404 ...
        }, ...
        'Rz_bar', { ...
            0.465,0.449,0.428,0.400,0.409,0.404,0.444,0.488,0.563,0.400,0.438,0.462,0.507,0.467,0.548,0.520,0.516,0.532,0.432,0.567,0.454,0.414,0.484 ...
        } ...
    );

    n = numel(db);

    % -----------------------------
    % Score candidates
    % Lower score = better match
    %
    % Uses normalized error in:
    %   - gross weight
    %   - span
    %   - length
    %
    % We weight geometry a bit more than weight.
    % -----------------------------
    scores = zeros(n,1);

    for i = 1:n
        wt_err = abs(db(i).GW_lb - MTOW_lb) / MTOW_lb;
        b_err  = abs(db(i).b_ft  - b_ft)    / b_ft;
        L_err  = abs(db(i).L_ft  - L_ft)    / L_ft;

        % weighted score
        scores(i) = 0.30*wt_err + 0.35*b_err + 0.35*L_err;
    end

    [~, idx] = min(scores);
    ref = db(idx);

    % -----------------------------
    % Dimensional radii of gyration
    % -----------------------------
    Rx_ft = ref.Rx_bar * b_ft;
    Ry_ft = ref.Ry_bar * L_ft;
    Rz_ft = ref.Rz_bar * b_ft;

    % -----------------------------
    % Inertia at MTOW
    % -----------------------------
    Ixx_MTOW = (MTOW_lb / g) * Rx_ft^2;
    Iyy_MTOW = (MTOW_lb / g) * Ry_ft^2;
    Izz_MTOW = (MTOW_lb / g) * Rz_ft^2;

    % -----------------------------
    % Inertia at empty weight
    % -----------------------------
    Ixx_empty = (Wempty_lb / g) * Rx_ft^2;
    Iyy_empty = (Wempty_lb / g) * Ry_ft^2;
    Izz_empty = (Wempty_lb / g) * Rz_ft^2;

    % -----------------------------
    % Output struct
    % -----------------------------
    inertia = struct();

    aircraft.inertia.reference_aircraft = ref.name;
    aircraft.inertia.reference_GW_lb = ref.GW_lb;
    aircraft.inertia.reference_b_ft = ref.b_ft;
    aircraft.inertia.reference_L_ft = ref.L_ft;

    aircraf.tinertia.Rx_bar = ref.Rx_bar;
    aircraft.inertia.Ry_bar = ref.Ry_bar;
    aircraft.inertia.Rz_bar = ref.Rz_bar;

    aircraft.inertia.Rx_ft = Rx_ft;
    aircraft.inertia.Ry_ft = Ry_ft;
    aircraft.inertia.Rz_ft = Rz_ft;

    aircraft.inertia.MTOW.Ixx_slug_ft2 = Ixx_MTOW;
    aircraft.inertia.MTOW.Iyy_slug_ft2 = Iyy_MTOW;
    aircraft.inertia.MTOW.Izz_slug_ft2 = Izz_MTOW;

    aircraft.inertia.Empty.Ixx_slug_ft2 = Ixx_empty;
    aircraft.inertia.Empty.Iyy_slug_ft2 = Iyy_empty;
    aircraft.inertia.Empty.Izz_slug_ft2 = Izz_empty;

    inertia.selection_score = scores(idx);
end