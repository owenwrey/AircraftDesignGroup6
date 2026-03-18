function [aircraft] = CgInertiaCalc(aircraft)
%CGINERTIACALC Calculates CG and products of inertia of aircraft Struct.
%   Detailed explanation will go here

%% main loop

%% calculate fuselage cg

aircraft.fuselage.cg.x = 1 + aircraft.fuselage.length/2;
aircraft.fuselage.cg.y = 0;
aircraft.fuselage.cg.z = 1 + aircraft.gear.mg.height + aircraft.fuselage.diameter/6;


%% calculate Engine cg's

aircraft.engine.cg.x = aircraft.fuselage.length - 15.192./2 + 1; % 1ft front, 1/2 of engine length
aircraft.engine.cg.y = 0;
aircraft.engine.cg.z = aircraft.fuselage.cg.z;

%% calculate avionics cg

aircraft.avionics.cg.x = 8;
aircraft.avionics.cg.y = 0;
aircraft.avionics.cg.z = aircraft.fuselage.cg.z + aircraft.fuselage.diameter/3;

%% calculate fuel tanks cg

combLength = aircraft.fuselage.length - 15.192./2 + 1 - 8;

aircraft.fuel1.cg.x = aircraft.cg.x-combLength/4;
aircraft.fuel1.cg.y = 0;
aircraft.fuel1.cg.z = aircraft.fuselage.cg.z;
aircraft.fuel1.weight = aircraft.weight.fuel/2;

aircraft.fuel2.cg.x = aircraft.cg.x+combLength/4;
aircraft.fuel2.cg.y = 0;
aircraft.fuel2.cg.z = aircraft.fuselage.cg.z;
aircraft.fuel2.weight = aircraft.weight.fuel/2;

%% calc wing cg



rt = 0; % root tip
rl = aircraft.wing.chord.root; % root tail
tt = -aircraft.wing.chord.tip/4 + (aircraft.wing.span/2)*tand(aircraft.wing.sweep); % tip tip
tl = aircraft.wing.chord.tip*(3/4) + (aircraft.wing.span/2)*tand(aircraft.wing.sweep); % tip tail

span = aircraft.wing.span/2;

aircraft.wing.poly = polyshape([rt rl tl tt],[0 0 span/2 span/2]);
[centX,~] = centroid(aircraft.wing.poly);


taper = aircraft.wing.taper_ratio;
ybar = (span/6)*((1 + 2*taper)/(1 + taper));

xglobcent = centX - aircraft.wing.chord.root/4 - ybar*tand(aircraft.wing.sweep) - aircraft.wing.MAC*(0.05) + aircraft.cg.x;

aircraft.wing.cg.x = xglobcent;
aircraft.wing.cg.y = 0;
aircraft.wing.cg.z = aircraft.fuselage.cg.z;

temp = aircraft.cg.x + aircraft.wing.MAC;


%% calc ordinance cg
aircraft.ordinance.cg.x = aircraft.cg.x - aircraft.wing.chord.root/4 - ybar*tand(aircraft.wing.sweep) - aircraft.wing.MAC*(0.05) + aircraft.wing.chord.root/4;
aircraft.ordinance.cg.y = 0;
aircraft.ordinance.cg.z = aircraft.wing.cg.z - 1;


%% calc Horiz stab cg



rt = 0; % root tip
rl = aircraft.ht.chord.root; % root tail
tt = -aircraft.ht.chord.tip/4 + (aircraft.ht.span/2)*tand(aircraft.ht.sweep); % tip tip
tl = aircraft.ht.chord.tip*(3/4) + (aircraft.ht.span/2)*tand(aircraft.ht.sweep); % tip tail

span = aircraft.ht.span/2;

aircraft.ht.poly = polyshape([rt rl tl tt],[0 0 span/2 span/2]);
[centX,~] = centroid(aircraft.ht.poly);


taper = aircraft.ht.TaperRatio;
ybar = (span/6)*((1 + 2*taper)/(1 + taper));

xglobcent = centX - aircraft.ht.chord.root/4 - ybar*tand(aircraft.ht.sweep) - aircraft.ht.MAC/2 + aircraft.cg.x + aircraft.ht.leverArm;

aircraft.ht.cg.x = xglobcent;
aircraft.ht.cg.y = 0;
aircraft.ht.cg.z = aircraft.fuselage.cg.z;


%% calc vert stab cg



rt = 0; % root tip
rl = aircraft.vt.chord.root_each; % root tail
tt = -aircraft.vt.chord.tip_each/4 + (aircraft.vt.span_each/2)*tand(aircraft.vt.sweep); % tip tip
tl = aircraft.vt.chord.tip_each*(3/4) + (aircraft.vt.span_each/2)*tand(aircraft.vt.sweep); % tip tail

span = aircraft.vt.span_each/2;

aircraft.vt.poly = polyshape([rt rl tl tt],[0 0 span/2 span/2]);
[centX,centY] = centroid(aircraft.vt.poly);


taper = aircraft.vt.TaperRatio;
ybar = (span/6)*((1 + 2*taper)/(1 + taper));

xglobcent = centX - aircraft.vt.chord.root_each/4 - ybar*tand(aircraft.vt.sweep) - aircraft.vt.MAC/2 + aircraft.cg.x + aircraft.vt.leverArm;



aircraft.vt.cg.x = xglobcent;
aircraft.vt.cg.y = 0;
aircraft.vt.cg.z = aircraft.fuselage.cg.z + centY*tand(45); % assume that vert stabs are at 45* angle

%% calc cockpit cg

aircraft.cockpit.cg.x = 8;
aircraft.cockpit.cg.y = 0;
aircraft.cockpit.cg.z = aircraft.fuselage.cg.z + aircraft.fuselage.diameter/3;

%% calc gear cg

aircraft.gear.mg.cg.x = aircraft.gear.mg.x;
aircraft.gear.mg.cg.y = 0;
aircraft.gear.mg.cg.z = aircraft.gear.mg.height/2 + 1;

aircraft.gear.ng.cg.x = aircraft.gear.ng.x;
aircraft.gear.ng.cg.y = 0;
aircraft.gear.ng.cg.z = aircraft.gear.mg.height/2 + 1;



%% calculate Center of Gravity
% place fields that shouldnt be iterated over here.
blacklist = ["flightcond","cg","constants","aero","enginesystems","weight","fuelSys","gear","TimeStepTable"];

% initialize variables
xsum = 0;
ysum = 0;
zsum = 0;
weightsum = 0;

field = fieldnames(aircraft)';

for i = 1:numel(field)
    
    if ~ismember(blacklist,field{i})

        
        % sum weights times locations
        xsum = aircraft.(field{i}).weight.*aircraft.(field{i}).cg.x + xsum;
        ysum = aircraft.(field{i}).weight.*aircraft.(field{i}).cg.y + ysum;
        zsum = aircraft.(field{i}).weight.*aircraft.(field{i}).cg.z + zsum;
        
        % sum weight
        weightsum = aircraft.(field{i}).weight + weightsum;
       
    end
end

%landing gear

    % NG sum weights times locations
        xsum = aircraft.gear.mg.weight.*aircraft.gear.mg.cg.x + xsum;
        ysum = aircraft.gear.mg.weight.*aircraft.gear.mg.cg.y + ysum;
        zsum = aircraft.gear.mg.weight.*aircraft.gear.mg.cg.z + zsum;
        
        % sum weight
        weightsum = aircraft.gear.mg.weight + weightsum;

   % MG sum weights times locations
        xsum = aircraft.gear.ng.weight.*aircraft.gear.ng.cg.x + xsum;
        ysum = aircraft.gear.ng.weight.*aircraft.gear.ng.cg.y + ysum;
        zsum = aircraft.gear.ng.weight.*aircraft.gear.ng.cg.z + zsum;
        
        % sum weight
        weightsum = aircraft.gear.ng.weight + weightsum;
        

% cg  & add to output
aircraft.cg.x = xsum./weightsum;
aircraft.cg.y = ysum./weightsum;
aircraft.cg.z = zsum./weightsum;



end