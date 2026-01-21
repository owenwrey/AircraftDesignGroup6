function comp = buildComponents(k_default)
% comp = buildComponents(k_default)

% wing
comp(1).name = 'wing';
comp(1).l_ft = 6.128;
comp(1).swet_ft2 = 1180.26;
comp(1).t_c = 0.045;
comp(1).x_c = 0.30;
comp(1).sweep_deg = 24;
comp(1).q = 1.0;
comp(1).k_ft = k_default;
comp(1).ff_scalar = []; % indicates wing-type FF computed

% fuselage
comp(2).name = 'fuselage';
comp(2).l_ft = 50;
comp(2).d_ft = 12;
comp(2).swet_ft2 = 1474.52;
f = comp(2).l_ft / comp(2).d_ft;
comp(2).ff_scalar = 0.9 + 5/f^1.5 + f/400;
comp(2).q = 1.0;
comp(2).k_ft = k_default;

% horizontal tail
comp(3).name = 'horizontal tail';
comp(3).l_ft = 5.785;
comp(3).swet_ft2 = 240.46;
comp(3).t_c = 0.04;
comp(3).x_c = 0.30;
comp(3).sweep_deg = 30;
comp(3).q = 1.03;
comp(3).k_ft = k_default;
comp(3).ff_scalar = [];

% vertical tail L
comp(4).name = 'vertical tail L';
comp(4).l_ft = 7.966;
comp(4).swet_ft2 = 116.10;
comp(4).t_c = 0.30;
comp(4).x_c = 0.30;
comp(4).sweep_deg = 35;
comp(4).q = 1.08;
comp(4).k_ft = k_default;
comp(4).ff_scalar = [];

% vertical tail R
comp(5) = comp(4);
comp(5).name = 'vertical tail R';

% strut
comp(6).name = 'strut';
comp(6).l_ft = 3;
comp(6).swet_ft2 = 5;
comp(6).t_c = 0.12;
comp(6).x_c = 0.30;
comp(6).sweep_deg = 5;
comp(6).q = 1.3;
comp(6).k_ft = k_default;
comp(6).ff_scalar = [];

% pylon
comp(7).name = 'pylon';
comp(7).l_ft = 4;
comp(7).swet_ft2 = 8;
comp(7).t_c = 0.12;
comp(7).x_c = 0.30;
comp(7).sweep_deg = 0;
comp(7).q = 1.4;
comp(7).k_ft = k_default;
comp(7).ff_scalar = [];

% nacelle
comp(8).name = 'nacelle';
comp(8).l_ft = 12;
comp(8).d_ft = 3.2;
comp(8).swet_ft2 = 120;
f = comp(8).l_ft / comp(8).d_ft;
comp(8).ff_scalar = 1 + 0.35/f;
comp(8).q = 1.5;
comp(8).k_ft = k_default;

end
