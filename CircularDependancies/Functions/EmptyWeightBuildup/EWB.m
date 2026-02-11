function weights = EWB(aircraft)

% initialize variables
weight.fuel = aircraft.constants.fuelWeight;
W0 = aircraft.constants.totalWeight;
weight.payload = aircraft.strike.weight;
weight.mg = aircraft.gear.mg.weight;
weight.ng = aircraft.gear.ng.weight;
weight.allElse = aircraft.avionics.weight;

cruiseDynamicPressure = aircraft.constants.cDP;
T2C = aircraft.wing.T2C;                    % thickness to chord ratio
N_lim = aircraft.constants.limitLoad;       % limit load factor (8g)
Nz = N_lim*1.5;                             % ultimate load factor

wing.area = aircraft.wing.area;
wing.AR = aircraft.wing.AR;
wing.QuarterChordSweep = aircraft.wing.sweep;
wing.TaperRatio = aircraft.wing.taper_ratio;
wing.span = aircraft.wing.span;
rootChord = aircraft.wing.rootChord;
tipChord = aircraft.wing.tipChord;
chordFrac = .2;             % control surface chord fraction

% calculate control surface areas
% ailerons: 20% local chord, going from 50% to 90% span
% flaps: 20% local chord, going from 0% to 50% span
chord50 = (rootChord-tipChord)*.5 + tipChord;       % chord @ 50% span
chord90 = (rootChord-tipChord)*.9 + tipChord;       % chord @ 90% span
chord0_20 = rootChord*chordFrac;                    % 20% chord @ root
chord50_20 = chord50*chordFrac;                     % 20% chord @ 50% span
chord90_20 = chord90*chordFrac;                     % 20% chord @ 90% span
span50 = wing.span*.5;                              % 50% wingspan
span40 = wing.span*.4;                              % 40% wingspan
flapArea = .5*(chord0_20 + chord50_20)*span50;      % flap area, ft^2
aileronArea = .5*(chord90_20 + chord50_20)*span40;  % aileron area, ft^2
S_csw = flapArea + aileronArea;                     % total control surface area


ht.area = aircraft.ht.Area;
ht.QuarterChordSweep = aircraft.ht.sweep;
ht.AspectRatio = aircraft.ht.AR;
ht.TaperRatio = aircraft.ht.TaperRatio;

vt.area = aircraft.vt.Area;
vt.QuarterChordSweep = aircraft.vt.sweep;
vt.AspectRatio = aircraft.vt.AR;
vt.TaperRatio = aircraft.vt.TaperRatio;


% calculate new weights
% from Raymer Fighter/Attack Weights, Sec 15.3.1
aircraft.wing.weight = .0103*(W0*Nz)^(.5) * (wing.area)^(.622) * (wing.AR)^(.785)...
                        * T2C * (1+wing.TaperRatio)^(.05) * ...
                        (cosd(wing.QuarterChordSweep))^(-1) * (S_csw)^(.04);




aircraft.ht.weight = .016*(ultimateLoadFactor*weight.total)^(.414) * cruiseDynamicPressure^(.168) * ht.area^(.896) * (100*thicknessToChordRatio/(cosd(ht.QuarterChordSweep)))^(-.12)...
    * (ht.AspectRatio/(cosd(ht.QuarterChordSweep))^2)^(.043) * ht.TaperRatio^(-.02);

aircraft.vt.weight = .073*(ultimateLoadFactor*weight.total)^(.376) * cruiseDynamicPressure^(.122) * vt.area^(.873) * (100*thicknessToChordRatio/(cosd(vt.QuarterChordSweep)))^(-.49)...
    * (vt.AspectRatio/(cosd(vt.QuarterChordSweep))^2)^(.357) * vt.TaperRatio^(.039);




prevWeight = weight.total;

weight.total = weight.fuel + weight.payload + weight.ng + weight.mg + weight.allElse + aircraft.wing.weight + aircraft.ht.weight + aircraft.vt.weight;

aircraft.constants.weight = weight.total;

weights = [prevWeight;aircraft.constants.weight];



end