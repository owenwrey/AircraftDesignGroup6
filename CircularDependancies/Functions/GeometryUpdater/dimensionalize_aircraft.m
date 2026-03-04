function aircraft = dimensionalize_aircraft(aircraft)
%DIMENSIONALIZE_AIRCRAFT  Compute aircraft dimensions from W0 and W/S.

%% Required checks 
req = @(s,f) isfield(s,f) && ~isempty(s.(f));

if ~isfield(aircraft,'constants') || ~req(aircraft.constants,'wingLoading')
    error('Need aircraft.constants.wingLoading.');
end
if ~isfield(aircraft,'weight') || ~req(aircraft.weight,'total')
    error('Need aircraft.weight.total.');
end
if ~isfield(aircraft,'fuselage') || ~req(aircraft.fuselage,'length') || ~req(aircraft.fuselage,'diameter')
    error('Need aircraft.fuselage.length and aircraft.fuselage.diameter.');
end
if ~isfield(aircraft,'wing') || ~req(aircraft.wing,'AR') || ~req(aircraft.wing,'taper_ratio')
    error('Need aircraft.wing.AR and aircraft.wing.taper_ratio.');
end
if ~isfield(aircraft,'ht') || ~req(aircraft.ht,'VolCoeff') || ~req(aircraft.ht,'AR') || ...
        ~req(aircraft.ht,'TaperRatio') || ~req(aircraft.ht,'leverArm_frac')
    error('Need aircraft.ht.VolCoeff, AR, TaperRatio, leverArm_frac.');
end
if ~isfield(aircraft,'vt') || ~req(aircraft.vt,'VolCoeff') || ~req(aircraft.vt,'AR') || ...
        ~req(aircraft.vt,'TaperRatio') || ~req(aircraft.vt,'leverArm_frac')
    error('Need aircraft.vt.VolCoeff, AR, TaperRatio, leverArm_frac.');
end

if ~isfield(aircraft.vt,'twinTail')
    aircraft.vt.twinTail = true;
end

%% Pull inputs
W0 = aircraft.weight.total;
WS = aircraft.constants.wingLoading;

Lf = aircraft.fuselage.length;
df = aircraft.fuselage.diameter;

ARw = aircraft.wing.AR;
TRw = aircraft.wing.taper_ratio;
T2C = aircraft.wing.T2C;

CHT = aircraft.ht.VolCoeff;
ARh = aircraft.ht.AR;
TRh = aircraft.ht.TaperRatio;

CVT = aircraft.vt.VolCoeff;
ARv = aircraft.vt.AR;
TRv = aircraft.vt.TaperRatio;

Lh = aircraft.ht.leverArm_frac * Lf;
Lv = aircraft.vt.leverArm_frac * Lf;

%% Fuselage 
aircraft.fuselage.volume = pi*(df/2)^2 * Lf;

%% Wing 
S_ref = W0 / WS;
b_w = sqrt(ARw * S_ref);
S_wet = 2*S_ref*(1 + 0.25 * T2C);

Croot_w = (2 * S_ref) / (b_w * (1 + TRw));
Ctip_w  = TRw * Croot_w;

MAC_w = ((2/3) * Croot_w * (1 + TRw + TRw^2)) / (1 + TRw);

aircraft.wing.Area = S_ref;
aircraft.wing.span = b_w;
aircraft.wing.swet = S_wet;
aircraft.wing.MAC  = MAC_w;
aircraft.wing.chord.root = Croot_w;
aircraft.wing.chord.tip  = Ctip_w;

%% Horizontal Tail 
S_HT = (CHT * MAC_w * S_ref) / Lh;
b_HT = sqrt(ARh * S_HT);


Croot_HT = (2 * S_HT) / (b_HT * (1 + TRh));
Ctip_HT  = TRh * Croot_HT;

MAC_HT = ((2/3) * Croot_HT * (1 + TRh + TRh^2)) / (1 + TRh);

aircraft.ht.leverArm = Lh;
aircraft.ht.Area     = S_HT;
aircraft.ht.span     = b_HT;
aircraft.ht.MAC      = MAC_HT;
aircraft.ht.chord.root = Croot_HT;
aircraft.ht.chord.tip  = Ctip_HT;

%% Vertical Tail 
S_VT_total = (CVT * b_w * S_ref) / Lv;   

if aircraft.vt.twinTail
    S_VT_each = S_VT_total / 2;
else
    S_VT_each = S_VT_total;
end

b_VT_each = sqrt(ARv * S_VT_each);

Croot_VT_each = (2 * S_VT_each) / (b_VT_each * (1 + TRv));
Ctip_VT_each  = TRv * Croot_VT_each;

MAC_VT_each = ((2/3) * Croot_VT_each * (1 + TRv + TRv^2)) / (1 + TRv);

aircraft.vt.leverArm  = Lv;
aircraft.vt.Area      = S_VT_total;   
aircraft.vt.Area_each = S_VT_each;
aircraft.vt.span_each = b_VT_each;
aircraft.vt.MAC_each  = MAC_VT_each;
aircraft.vt.chord.root_each = Croot_VT_each;
aircraft.vt.chord.tip_each  = Ctip_VT_each;

end
