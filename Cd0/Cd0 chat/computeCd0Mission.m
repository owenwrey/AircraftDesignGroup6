function Cd0_total = computeCd0Mission(alt_ft, M_vec, comp, params)
% Cd0_total = computeCd0Mission(alt_ft, M_vec, comp, params)

n = numel(alt_ft);
Cd0_total = zeros(n,1);

for j = 1:n
    M = M_vec(j);

    if M < 1
        Cd0_total(j) = cd0Subsonic(alt_ft(j), M, comp, params);
    else
        Cd0_total(j) = cd0Supersonic(alt_ft(j), M, comp, params);
    end
end

end
