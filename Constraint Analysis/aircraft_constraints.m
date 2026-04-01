function [c, ceq] = aircraft_constraints(x, W2S_plot, T2W_Ceiling, T2W_MaxSpeed, WS_InstTurn)
    % x(1) is Wing Loading (W/S)
    % x(2) is Thrust-to-Weight (T/W)
    ceq = []; 

    % 1. Interpolate the curves to find the requirement at the current x(1)
    % We use 'linear' interpolation and 'extrap' to ensure the GA doesn't crash
    req_Ceiling = interp1(W2S_plot, T2W_Ceiling, x(1), 'linear', 'extrap');
    req_MaxSpeed = interp1(W2S_plot, T2W_MaxSpeed, x(1), 'linear', 'extrap');

    % 2. Define Inequalities (c <= 0 is the "Safe" zone)
    % We want T/W (x2) to be GREATER than the requirements
    c(1) = req_Ceiling - x(2); 
    c(2) = req_MaxSpeed - x(2);
    
    % 3. Vertical Limit (Instantaneous Turn)
    % Wing loading x(1) must be LESS than the limit
    c(3) = x(1) - WS_InstTurn;
end