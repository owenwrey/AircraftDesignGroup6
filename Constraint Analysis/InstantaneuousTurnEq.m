alt = 20000 * 0.3048;   
n = 7;
v = 295; 
Cl = 1.2;
beta = 0.9

[T_ISA,~,~,rho_ISA] = atmosisa(alt,"extended","on");
T = T_ISA + DeltaT;
rho = rho_ISA .* (T_ISA ./ T);

W2S_InstTurn = verticalConstraintAnalysis(n,beta,Cl,k,v,rho);

q = .5 .* rho .* v.^2;
WingLoading = q.*Cl./n./beta./k.^2
