function WingLoading = verticalConstraintAnalysis(n,beta,Cl,k,v,rho)
% Description: calculates wing loading for a given condition.
%
% n = load factor
% beta = weight lapse
% Cl = max lift coeficient
% k = safety factor
% v = velocity {m/s}
% rho = air density {kg/m^3}

q = .5 .* rho .* v.^2;
WingLoading = q.*Cl./n./beta./k.^2;
end