function [T, a, p, rho] = atmosAtFt(alt_ft)
% [T, a, p, rho] = atmosAtFt(alt_ft)

alt_m = alt_ft * 0.3048;
[T, a, p, rho] = atmosisa(alt_m);

end
