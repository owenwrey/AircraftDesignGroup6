function [a_mps, rho_slugft3] = atmosAtAlt(alt_ft, p)
%ATMOSATALT Wrapper around atmosisa; returns speed of sound (m/s) and rho (slug/ft^3).
[~, a_mps, ~, rho_kgm3] = atmosisa(alt_ft * p.ft2m);
rho_slugft3 = rho_kgm3 * 0.00194032;
end
