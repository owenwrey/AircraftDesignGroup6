function geom = compute_geometry(aircraft)
%COMPUTE_GEOMETRY  Computes geometry 

%% pull syntax 
W0 = aircraft.constants.W0;
WS = aircraft.constants.WS;

Lf = aircraft.fuselage.L;
df = aircraft.fuselage.d;

ARw = aircraft.wing.AR;
TRw = aircraft.wing.taper_ratio;

CHT = aircraft.ht.C;
ARh = aircraft.ht.AR;
TRh = aircraft.ht.taper_ratio;
Lh  = aircraft.ht.L_percent * Lf;

CVT = aircraft.vt.C;
ARv = aircraft.vt.AR;
TRv = aircraft.vt.taper_ratio;
Lv  = aircraft.vt.L_percent * Lf;

isTwin = isfield(aircraft.vt,'twin') && aircraft.vt.twin;

%% Fuselage 
geom.fuselage.L = Lf;
geom.fuselage.d = df;
geom.fuselage.volume = pi*(df/2)^2 * Lf;

%% Main Wing 
geom.wing.area = W0/WS;
geom.wing.span = sqrt(ARw * geom.wing.area);

geom.wing.croot = (2*geom.wing.area) / (geom.wing.span * (1 + TRw));
geom.wing.ctip  = TRw * geom.wing.croot;

geom.wing.MAC = (2/3) * geom.wing.croot * (1 + TRw + TRw^2) / (1 + TRw);
geom.wing.Y   = (geom.wing.span/6) * ((1 + 2*TRw) / (1 + TRw));

%% Horizontal Tail 
geom.ht.L = Lh;

geom.ht.area = (CHT * geom.wing.MAC * geom.wing.area) / geom.ht.L;
geom.ht.span = sqrt(ARh * geom.ht.area);

geom.ht.croot = (2*geom.ht.area) / (geom.ht.span * (1 + TRh));
geom.ht.ctip  = TRh * geom.ht.croot;

geom.ht.MAC = (2/3) * geom.ht.croot * (1 + TRh + TRh^2) / (1 + TRh);
geom.ht.Y   = (geom.ht.span/6) * ((1 + 2*TRh) / (1 + TRh));

%% Vertical Tail (twin tail) 
geom.vt.L = Lv;

geom.vt.area_total = (CVT * geom.wing.span * geom.wing.area) / geom.vt.L; % corrected form

if isTwin
    geom.vt.area_each = geom.vt.area_total / 2;
else
    geom.vt.area_each = geom.vt.area_total;
end

geom.vt.span_each = sqrt(ARv * geom.vt.area_each);

geom.vt.croot_each = (2*geom.vt.area_each) / (geom.vt.span_each * (1 + TRv));
geom.vt.ctip_each  = TRv * geom.vt.croot_each;

geom.vt.MAC_each = (2/3) * geom.vt.croot_each * (1 + TRv + TRv^2) / (1 + TRv);
geom.vt.Y_each   = (geom.vt.span_each/6) * ((1 + 2*TRv) / (1 + TRv));

end