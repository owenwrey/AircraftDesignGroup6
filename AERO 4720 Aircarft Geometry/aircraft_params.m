function aircraft = aircraft_params()

%% aircraft.constants 
aircraft.constants.W0 = 76965;   % [lb]
aircraft.constants.WS = 125;     % wing loading [lb/ft^2]

%% aircraft.fuselage
aircraft.fuselage.L = 50;        % [ft]
aircraft.fuselage.d = 6.5;       % [ft]

%% aircraft.wing 
aircraft.wing.AR = 4;
aircraft.wing.taper_ratio = 0.25;   
aircraft.wing.MAC  = NaN;            % output later, placeholder ok
aircraft.wing.span = NaN;            % output later
aircraft.wing.sweep = NaN;           % placeholder
aircraft.wing.weight = NaN;          % placeholder

%% aircraft.ht 
aircraft.ht.C = 0.4;                 % C_HT
aircraft.ht.AR = 4;
aircraft.ht.taper_ratio = 0.4;
aircraft.ht.L_percent = 0.525;       % fraction of fuselage length
aircraft.ht.CD0 = NaN;               % placeholder

% aircraft.vt 
aircraft.vt.C = 0.09;                % C_VT
aircraft.vt.AR = 1;
aircraft.vt.taper_ratio = 0.3;
aircraft.vt.L_percent = 0.475;
aircraft.vt.twin = true;

end