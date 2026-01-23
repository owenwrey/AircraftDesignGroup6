function Thrust_2_Weight = Master_Eqn(CD_0,k1,k2,CD_R,n,v,rho,ddt_h,ddt_v,alpha,beta,W2S)
% Description: calculates thrust to weight using the master equation.
%
% Cd0 = zero-lift drag coeficient
% K1 = drag polar constant 1
% K2 = drag polar constant 2
% CdR = resultant drag coeficient
% n = load factor
% v = velocity {m/s}
% rho = density {kg/m^3}
% ddt_h = climb speed {m/s}
% ddt_v = acceleration {m/s^2}
% alpha = thrust lapse
% beta = weight lapse
% W2S = wingloading {N/m^2}



g = 9.81; % {m/s^2}
q = .5 .*rho .* v.^2; % {Pa}

Thrust_2_Weight = beta ./ alpha .* (q./beta./W2S.*(CD_0 + k1.*(n.*beta.*W2S./q) + k2.*(n.*beta.*W2S./q).^2 + CD_R) + ddt_h/v + ddt_v/g);
end