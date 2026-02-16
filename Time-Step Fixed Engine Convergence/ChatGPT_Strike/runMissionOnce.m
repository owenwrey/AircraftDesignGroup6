function [Tbl, FuelReq, W0_next] = runMissionOnce(W0, FuelReq, SegNames, npts)
%RUNMISSIONONCE Runs a single mission evaluation given W0 and FuelReq guess.

p = initParams();
p.W0 = W0;
p.FuelReq = FuelReq;

[Tbl, segIdx] = initMissionTable(SegNames, npts);

% Segment 1 initializes Weight/FuelRem.
Tbl = seg01_SWT(Tbl, segIdx(1,1), segIdx(1,2), p);
Tbl = seg02_TKO(Tbl, segIdx(2,1), segIdx(2,2), p);
Tbl = seg03_CLIMB1(Tbl, segIdx(3,1), segIdx(3,2), p);
Tbl = seg04_CRUISE_OBD(Tbl, segIdx(4,1), segIdx(4,2), p);
Tbl = seg05_DESC1(Tbl, segIdx(5,1), segIdx(5,2), p);
Tbl = seg06_COMBAT(Tbl, segIdx(6,1), segIdx(6,2), p);
Tbl = seg07_WP_FIRE(Tbl, segIdx(7,1), segIdx(7,2), p);
Tbl = seg08_CLIMB2(Tbl, segIdx(8,1), segIdx(8,2), p);
Tbl = seg09_CRUISE_IBD(Tbl, segIdx(9,1), segIdx(9,2), p);
Tbl = seg10_DESC2(Tbl, segIdx(10,1), segIdx(10,2), p);
Tbl = seg11_LTR2(Tbl, segIdx(11,1), segIdx(11,2), p);
Tbl = seg12_LTS(Tbl, segIdx(12,1), segIdx(12,2), p);

% Post-processing (same intent as original)
Tbl.EnHt = Tbl.Alt + (Tbl.KTAS*p.NM2ft).^2/(2*32.17);

% Weight iteration block from original script
A = 2.34;
C = -0.13;
EWF = A * W0^C;

W_crew    = 300;
W_payload = p.payloadDrop_lb;

OEW = EWF * W0;

FuelAllow = 0.06 * Tbl.FuelBurn(segIdx(12,2)); % 6% fuel allowance
FuelReq   = FuelAllow + Tbl.FuelBurn(segIdx(12,2));

FuelAvail = W0 - OEW - W_crew - W_payload;

W0_next = FuelReq + OEW + W_crew + W_payload;

% (FuelAvail and FuelExcess are available to caller if desired)
end
