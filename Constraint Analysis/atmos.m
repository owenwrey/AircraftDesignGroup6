function [T,P,rho] = atmos(alt,del_T)
%[T,a,P,rho] = atmosisa(alt);
if nargin == 1
    del_T = 0;
end
            
T = 15.04  - 0.00649.*alt;
P = 101.29 .* ((T + 273.1)./288.08).^(5.256);
T = T + del_T;
rho = P ./ (.2869.*(T+273.1));
end