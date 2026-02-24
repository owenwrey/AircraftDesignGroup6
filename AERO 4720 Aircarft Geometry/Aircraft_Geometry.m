clc; clear;

aircraft = aircraft_params();          
geom    = compute_geometry(aircraft);  

aircraft.wing.MAC  = geom.wing.MAC;
aircraft.wing.span = geom.wing.span;

%% PRINT RESULTS
fprintf("\nFuselage\n");
fprintf("Fuselage Volume (ft^3): %.2f\n", geom.fuselage.volume);

fprintf("\nMain Wing\n");
fprintf("Wing Area (ft^2): %.2f\n", geom.wing.area);
fprintf("Wingspan (ft): %.2f\n", geom.wing.span);
fprintf("Root Chord (ft): %.2f\n", geom.wing.croot);
fprintf("Tip Chord  (ft): %.2f\n", geom.wing.ctip);
fprintf("MAC (ft): %.2f\n", geom.wing.MAC);
fprintf("Y (ft): %.2f\n", geom.wing.Y);

fprintf("\nHorizontal Tail\n");
fprintf("Area (ft^2): %.2f\n", geom.ht.area);
fprintf("Span (ft): %.2f\n", geom.ht.span);
fprintf("Root Chord (ft): %.2f\n", geom.ht.croot);
fprintf("Tip Chord  (ft): %.2f\n", geom.ht.ctip);
fprintf("MAC (ft): %.2f\n", geom.ht.MAC);
fprintf("Y (ft): %.2f\n", geom.ht.Y);

fprintf("\nVertical Tail\n");
fprintf("Total Area (ft^2): %.2f\n", geom.vt.area_total);
fprintf("Each Tail Area (ft^2): %.2f\n", geom.vt.area_each);
fprintf("Each Tail Span (ft): %.2f\n", geom.vt.span_each);
fprintf("Root Chord Each (ft): %.2f\n", geom.vt.croot_each);
fprintf("Tip Chord  Each (ft): %.2f\n", geom.vt.ctip_each);
fprintf("MAC Each (ft): %.2f\n", geom.vt.MAC_each);
fprintf("Y Each (ft): %.2f\n", geom.vt.Y_each);