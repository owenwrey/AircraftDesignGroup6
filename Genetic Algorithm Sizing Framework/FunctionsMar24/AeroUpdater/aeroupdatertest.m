clc;clear;close all

k_default = 0.50e-5
% altitude range 
alt = 0:500:25000;
% mach range
m = 0.1:0.1:2;

% wing 
aircraft.wing.name = 'wing';
aircraft.wing.l    = 6.128;        % ft (updated MAC from board)
aircraft.wing.swet = 1180.26;      % ft^2 wetted area from board geometry
aircraft.wing.t_c  = 0.055;        % thickness ratio t/c 
aircraft.wing.x_c  = 0.30;         % x/c location of max thickness
aircraft.wing.sweep = 30;    % degrees  (F/A-18 style LE sweep)
aircraft.wing.sweep = deg2rad(aircraft.wing.sweep);
aircraft.wing.ff = (1 + (0.6/aircraft.wing.x_c)*aircraft.wing.t_c + 100*(aircraft.wing.t_c)^4) .* ...
              (1.34*m.^0.18 .* cos(aircraft.wing.sweep).^0.28);
aircraft.wing.q = 1.0;         % pg 425, the f18 has a midwing, so q = 1
aircraft.wing.k = k_default;

% fuselage 
aircraft.fuselage.name = 'fuselage';
aircraft.fuselage.l    = 48;       % ft fuselage length (board)
aircraft.fuselage.d    = 12;       % ft max diameter (board)
aircraft.fuselage.swet = 1474.52;      % ft^2 computed wetted area from board dims
aircraft.fuselage.f    = aircraft.fuselage.l / aircraft.fuselage.d; 
aircraft.fuselage.ff   = 0.9 + 5/aircraft.fuselage.f^1.5 + aircraft.fuselage.f/400;
aircraft.fuselage.q    = 1.0;      % pg 429, the fuselage usually has negligible q
aircraft.fuselage.k    = k_default;


% horizontal tail
aircraft.ht.name = 'horizontal tail';
aircraft.ht.l    = 5.785;       % ft MAC from board
aircraft.ht.swet = 240.46;      % ft^2 wetted area from board geometry
aircraft.ht.t_c  = 0.05;
aircraft.ht.x_c  = 0.30;
aircraft.ht.sweep = 30;    % realistic HT sweep for fighter
aircraft.ht.sweep = deg2rad(aircraft.ht.sweep);
aircraft.ht.ff = (1 + (0.6/aircraft.ht.x_c)*aircraft.ht.t_c + 100*(aircraft.ht.t_c)^4) .* ...
              (1.34*m.^0.18 .* cos(aircraft.ht.sweep).^0.28);
aircraft.ht.q = 1.03;        % pg 429
aircraft.ht.k = k_default;

% vertical tail (left)

aircraft.vt.left.name = 'vertical tail L';
aircraft.vt.left.l    = 7.966;      % ft MAC (from board VT geometry)
aircraft.vt.left.swet = 116.10;     % ft^2 wetted area per tail (from board)
aircraft.vt.left.t_c  = 0.045;
aircraft.vt.left.x_c  = 0.30;
aircraft.vt.left.sweep = 35;
aircraft.vt.left.sweep= deg2rad(aircraft.vt.left.sweep);
aircraft.vt.left.ff = (1 + (0.6/aircraft.vt.left.x_c)*aircraft.vt.left.t_c + 100*(aircraft.vt.left.t_c)^4) .* ...
              (1.34*m.^0.18 .* cos(aircraft.vt.left.sweep).^0.28);
aircraft.vt.left.q = 1.08;         % pg 429
aircraft.vt.left.k = k_default;

% vertical tail (right)
aircraft.vt.right.name = 'vertical tail R';
aircraft.vt.right.l    = 7.966;      % ft MAC
aircraft.vt.right.swet = 116.10;     % ft^2 wetted area
aircraft.vt.right.t_c  = 0.045;
aircraft.vt.right.x_c  = 0.30;
aircraft.vt.right.sweep = 35;
aircraft.vt.right.sweep = deg2rad(aircraft.vt.right.sweep);
aircraft.vt.right.ff = (1 + (0.6/aircraft.vt.right.x_c)*aircraft.vt.right.t_c + 100*(aircraft.vt.right.t_c)^4) .* ...
              (1.34*m.^0.18 .* cos(aircraft.vt.right.sweep).^0.28);
aircraft.vt.right.q = 1.08;         % pg 429
aircraft.vt.right.k = k_default;

% strut
aircraft.strut.name = 'strut';
aircraft.strut.l    = 3;        % ft
aircraft.strut.swet = 5;        % ft^2
aircraft.strut.t_c  = 0.12;
aircraft.strut.x_c  = 0.30;
aircraft.strut.sweep = 5;
aircraft.strut.sweep = deg2rad(aircraft.strut.sweep);
aircraft.strut.ff = (1 + (0.6/aircraft.strut.x_c)*aircraft.strut.t_c + 100*(aircraft.strut.t_c)^4) .* ...
              (1.34*m.^0.18 .* cos(aircraft.strut.sweep).^0.28);
aircraft.strut.q = 1.3;         % pg 429, wing strut has < drag than pylon
aircraft.strut.k = k_default;

% pylon 
aircraft.pylon.name = 'pylon';
aircraft.pylon.l    = 4;        % ft
aircraft.pylon.swet = 8;        % ft^2
aircraft.pylon.t_c  = 0.12;
aircraft.pylon.x_c  = 0.30;
aircraft.pylon.sweep = 0;
aircraft.pylon.sweep = deg2rad(aircraft.pylon.sweep);
aircraft.pylon.ff = (1 + (0.6/aircraft.pylon.x_c)*aircraft.pylon.t_c + 100*(aircraft.pylon.t_c)^4) .* ...
              (1.34*m.^0.18 .* cos(aircraft.pylon.sweep).^0.28);
aircraft.pylon.q = 1.4;         % pg 429
aircraft.pylon.k = k_default;

% nacelle
aircraft.nacelle.name = 'nacelle';
aircraft.nacelle.l    = 12;       % ft
aircraft.nacelle.d    = 3.2;      % ft
aircraft.nacelle.swet = 120;      % ft^2
aircraft.nacelle.f    = aircraft.nacelle.l / aircraft.nacelle.d;
aircraft.nacelle.ff   = 1 + 0.35/aircraft.nacelle.f;
aircraft.nacelle.q    = 1.5;      % pg 425 - nacelle mounted directly on wing
aircraft.nacelle.k    = k_default;

aircraft = aeroupdater(aircraft);
