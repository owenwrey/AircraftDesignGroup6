function [Tbl, segIdx] = initMissionTable(SegNames, npts)
%INITMISSIONTABLE Pre-allocates the mission table and segment index ranges.
%
% Outputs:
%   Tbl    table with preallocated columns
%   segIdx Nx2 array with [istart iend] row for each segment

nSeg = numel(SegNames);
npts_sum = sum(npts);

segIdx = zeros(nSeg,2);
ctr = 1;
for s = 1:nSeg
    segIdx(s,1) = ctr;
    segIdx(s,2) = ctr + npts(s) - 1;
    ctr = ctr + npts(s);
end

Tbl = table();
Tbl.SegNum = zeros(npts_sum,1);
Tbl.Seg    = repmat({' '},npts_sum,1);

for s = 1:nSeg
    Tbl.Seg(segIdx(s,1):segIdx(s,2)) = SegNames(s);
    Tbl.SegNum(segIdx(s,1):segIdx(s,2)) = s;
end

% Kinematics / state
Tbl.Time  = zeros(npts_sum,1); % min
Tbl.dTime = zeros(npts_sum,1); % min
Tbl.Alt   = zeros(npts_sum,1); % ft
Tbl.rho   = zeros(npts_sum,1); % slug/ft^3

Tbl.KEAS  = zeros(npts_sum,1); % kt
Tbl.KTAS  = zeros(npts_sum,1); % kt
Tbl.MACH  = zeros(npts_sum,1); % -
Tbl.GS    = zeros(npts_sum,1); % kt
Tbl.FPA   = zeros(npts_sum,1); % deg
Tbl.EnHt  = zeros(npts_sum,1); % ft

Tbl.Dist  = zeros(npts_sum,1); % NM
Tbl.dDist = zeros(npts_sum,1); % NM
Tbl.dhdt  = zeros(npts_sum,1); % ft/min
Tbl.dVdt  = zeros(npts_sum,1); % ft/s^2 (as used in original)

% Aero / propulsion
Tbl.Weight = zeros(npts_sum,1); % lb
Tbl.WtDrop = zeros(npts_sum,1); % lb
Tbl.WtFrac = zeros(npts_sum,1); % -

Tbl.CL   = zeros(npts_sum,1);
Tbl.CD0  = zeros(npts_sum,1);
Tbl.K1   = zeros(npts_sum,1);
Tbl.K2   = zeros(npts_sum,1);
Tbl.CDR  = zeros(npts_sum,1);
Tbl.CD   = zeros(npts_sum,1);
Tbl.L_D  = zeros(npts_sum,1);

Tbl.Drag   = zeros(npts_sum,1); % lbf
Tbl.Thrust = zeros(npts_sum,1); % lbf
Tbl.TLapse = zeros(npts_sum,1); % -
Tbl.Ps     = zeros(npts_sum,1); % ft/min
Tbl.THROT  = zeros(npts_sum,1); % -

% Fuel accounting
Tbl.FF       = zeros(npts_sum,1); % lb/hr
Tbl.dFuel    = zeros(npts_sum,1); % lb
Tbl.FuelRem  = zeros(npts_sum,1); % lb
Tbl.FuelBurn = zeros(npts_sum,1); % lb
end
