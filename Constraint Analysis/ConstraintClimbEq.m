alt = 0;
n = 1;
M = 0.8;

[T_ISA,a_ISA,~,rho_ISA] = atmosisa(alt,"extended","on");
T = T_ISA + DeltaT;
a = sqrt(gamma*R*T);
rho = rho_ISA .* (T_ISA ./ T);

v = M.*a;
ddt_h = 28000*0.3048/60;
ddt_v = 0;
alpha = 1.11*(rho/rho_SL)-0.11;
beta = 1.0;

g = 9.81; % {m/s^2}
q = .5 .*rho .* v.^2;

Thrust_2_Weight = 1 / 0.9382 .* (45393/1/W2S.*(0.02 + + 0.05.*(1.*1*W2S./45393).^2) + 142.24/280)
