%GA Equations

alt = 30000*0.3048;
n = 1;
M = 1.6;

[T_ISA,a_ISA,~,rho_ISA] = atmosisa(alt,"extended","on");
T = T_ISA + DeltaT;
a = sqrt(gamma*R*T);
rho = rho_ISA .* (T_ISA ./ T);

v = M.*a;
ddt_h = 0;
ddt_v = 0;
alpha = (1.11*(rho/rho_SL)-0.11)*bump;
beta = 0.85;

g = 9.81; % {m/s^2}
q = .5 .*rho .* v.^2;

Thrust_2_Weight = 0.85 / 0.3733 .* (53920./0.85./W2S.*(0.02 + 0.05.*(1.*0.85.*W2S./53920).^2))
