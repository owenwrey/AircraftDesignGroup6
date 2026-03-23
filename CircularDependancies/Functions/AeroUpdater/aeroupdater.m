function aircraft = aeroupdater(Aircraft)
%% AEROUPDATER
% Builds C_D0 as a function of:
%   1) altitude, h [ft]
%   2) Mach number, M [-]
%
% Stores results in:
%   Aircraft.constants.alt_grid_ft
%   Aircraft.constants.mach_grid
%   Aircraft.constants.cd0_total_grid
%   Aircraft.constants.cd0_sum_grid
%   Aircraft.constants.cd0_misc_grid
%   Aircraft.constants.cd0_wave_grid
%   Aircraft.constants.cd0_lp_grid
%
% Also creates:
%   Aircraft.aero.cd0_interp
%   Aircraft.aero.cd0
%
% So later Jordan can call:
%   cd0_val = Aircraft.aero.cd0(alt_ft, mach)

%% -------------------------
% USER INPUTS
% --------------------------
alt_ft = 0:1000:50000;      % altitude grid [ft]
mach   = 0.1:0.1:2.0;       % Mach grid [-]

%% -------------------------
% CONSTANTS
% --------------------------
k_default = 0.50e-5;              % equivalent roughness height [ft]
sref      = Aircraft.wing.Area;  % reference wing area [ft^2]

%% -------------------------
% COMPONENT BUILDUP
% --------------------------
% Wing
comp(1).name = 'wing';
comp(1).l    = Aircraft.wing.MAC;          % [ft]
comp(1).swet = Aircraft.wing.swet;         % [ft^2]
comp(1).t_c  = Aircraft.wing.T2C;          % [-]
comp(1).x_c  = Aircraft.wing.x_c;          % [-]
comp(1).sweep_angle_deg = Aircraft.wing.sweep;
comp(1).sweep_angle     = deg2rad(comp(1).sweep_angle_deg);
comp(1).ff = (1 + (0.6/comp(1).x_c)*comp(1).t_c + 100*(comp(1).t_c)^4) .* ...
             (1.34*mach.^0.18 .* cos(comp(1).sweep_angle).^0.28);
comp(1).q = 1.00;
comp(1).k = k_default;

% Fuselage
comp(2).name = 'fuselage';
comp(2).l    = Aircraft.fuselage.length;   % [ft]
comp(2).d    = Aircraft.fuselage.diameter; % [ft]
comp(2).swet = Aircraft.fuselage.swet;     % [ft^2]
comp(2).f    = comp(2).l / comp(2).d;      % fineness ratio [-]
comp(2).ff   = 0.9 + 5/(comp(2).f^1.5) + comp(2).f/400;
comp(2).q    = 1.00;
comp(2).k    = k_default;

% Horizontal tail
comp(3).name = 'horizontal tail';
comp(3).l    = Aircraft.ht.MAC;            % [ft]
comp(3).swet = Aircraft.ht.swet;           % [ft^2]
comp(3).t_c  = Aircraft.ht.T2C;            % [-]
comp(3).x_c  = Aircraft.ht.x_c;            % [-]
comp(3).sweep_angle_deg = Aircraft.ht.sweep;
comp(3).sweep_angle     = deg2rad(comp(3).sweep_angle_deg);
comp(3).ff = (1 + (0.6/comp(3).x_c)*comp(3).t_c + 100*(comp(3).t_c)^4) .* ...
             (1.34*mach.^0.18 .* cos(comp(3).sweep_angle).^0.28);
comp(3).q = 1.03;
comp(3).k = k_default;

% Vertical tail LEFT
comp(4).name = 'vertical tail L';
comp(4).l    = Aircraft.vt.MAC;            % [ft]
comp(4).swet = Aircraft.vt.swet;           % [ft^2]
comp(4).t_c  = Aircraft.vt.T2C;            % [-]
comp(4).x_c  = Aircraft.vt.x_c;            % [-]
comp(4).sweep_angle_deg = Aircraft.vt.sweep;
comp(4).sweep_angle     = deg2rad(comp(4).sweep_angle_deg);
comp(4).ff = (1 + (0.6/comp(4).x_c)*comp(4).t_c + 100*(comp(4).t_c)^4) .* ...
             (1.34*mach.^0.18 .* cos(comp(4).sweep_angle).^0.28);
comp(4).q = 1.08;
comp(4).k = k_default;

% Vertical tail RIGHT
comp(5).name = 'vertical tail R';
comp(5).l    = comp(4).l;
comp(5).swet = comp(4).swet;
comp(5).t_c  = comp(4).t_c;
comp(5).x_c  = comp(4).x_c;
comp(5).sweep_angle_deg = comp(4).sweep_angle_deg;
comp(5).sweep_angle     = comp(4).sweep_angle;
comp(5).ff   = comp(4).ff;
comp(5).q    = comp(4).q;
comp(5).k    = comp(4).k;

%% -------------------------
% PREALLOCATE STORAGE GRIDS
% --------------------------
nAlt  = length(alt_ft);
nMach = length(mach);

cd0_sum_grid   = zeros(nAlt, nMach);   % component buildup only
cd0_misc_grid  = zeros(nAlt, nMach);   % stores / misc drag
cd0_wave_grid  = zeros(nAlt, nMach);   % wave drag
cd0_lp_grid    = zeros(nAlt, nMach);   % leakage + protuberance
cd0_total_grid = zeros(nAlt, nMach);   % total C_D0

%% -------------------------
% MAIN LOOPS
% --------------------------
for j = 1:nAlt

    h_ft = alt_ft(j);          % [ft]
    h_m  = h_ft * 0.3048;      % [m]

    % Standard atmosphere in SI
    [T, a, ~, rho] = atmosisa(h_m);
    mu = sutherland_mu(T);     % [kg/(m*s)]

    for k = 1:nMach

        M = mach(k);           % [-]
        V = M * a;             % [m/s]

        cd0_sum  = 0;
        CD0_misc = 0;
        CD0_wave = 0;

        %% -------------------------
        % MISC / STORES DRAG
        % --------------------------
        A_base = [0.136 0.136 1.068];   % [ft^2] [AIM-120, AIM-9X, MK-83]

        if M > 1
            d_q = (0.064 + 0.042*(M - 3.84)^2) .* A_base;   % [ft^2]
        else
            d_q = (0.139 + 0.419*(M - 0.161)^2) .* A_base;  % [ft^2]
        end

        Cd_misc_air_to_air = (6*d_q(1) + 2*d_q(2)) / sref; %#ok<NASGU>
        Cd_misc_strike     = (4*d_q(3) + 2*d_q(2)) / sref;

        % Choose the loadout you want
        CD0_misc = Cd_misc_strike;

        %% -------------------------
        % COMPONENT FRICTION + FORM / INTERFERENCE BUILDUP
        % --------------------------
        for i = 1:length(comp)

            l_ft = comp(i).l;      % [ft]
            l_m  = l_ft * 0.3048;  % [m]
            swet = comp(i).swet;   % [ft^2]
            k_c  = comp(i).k;      % [ft]

            % Reynolds number
            Re = rho * V * l_m / mu;

            if M > 1
                % Supersonic: remove subsonic FF/Q corrections
                FF = 1.0;
                Q  = 1.0;

                Re_cut = 44.62 * (l_ft/k_c)^1.053 * M^1.16;
                Re_eff = min(Re, Re_cut);

                Cf = 0.455 / ( (log10(Re_eff))^2.58 * (1 + 0.144*M^2)^0.65 );
                cd0_i = (Cf * FF * Q * swet) / sref;

            else
                q_c = comp(i).q;

                if isscalar(comp(i).ff)
                    ff_i = comp(i).ff;
                else
                    ff_i = comp(i).ff(k);
                end

                Re_cut = 38.21 * (l_ft/k_c)^1.053;
                Re_eff = min(Re, Re_cut);

                Cf = 0.455 / ( (log10(Re_eff))^2.58 * (1 + 0.144*M^2)^0.65 );
                cd0_i = (Cf * ff_i * q_c * swet) / sref;
            end

            cd0_sum = cd0_sum + cd0_i;
        end

        %% -------------------------
        % WAVE DRAG
        % --------------------------
        if M > 1
            d_fus = comp(2).d;               % [ft]
            ell   = comp(2).l * 0.50;        % [ft]
            Amax  = pi * (d_fus/2)^2;        % [ft^2]
            Ewd   = 2.0;                     % empirical correction factor

            if M < 1.2
                CD0_wave = ((9*pi/2) * (Amax/ell)^2 * Ewd) / sref;
            else
                Dq_SH = (9*pi/2) * (Amax/ell)^2;   % [ft^2]
                Lambda_LE_deg = comp(1).sweep_angle_deg;

                corr = 1 - 0.2*(M - 1.2)^0.57 * (1 - (pi*Lambda_LE_deg^0.77)/100);
                Dq_wave  = Ewd * corr * Dq_SH;     % [ft^2]
                CD0_wave = Dq_wave / sref;
            end
        else
            CD0_wave = 0;
        end

        %% -------------------------
        % LEAKAGE + PROTUBERANCE DRAG
        % --------------------------
        cd0_lp = 0.12 * (cd0_sum + CD0_misc + CD0_wave);

        %% -------------------------
        % TOTAL C_D0
        % --------------------------
        cd0_total = cd0_sum + CD0_misc + CD0_wave + cd0_lp;

        %% -------------------------
        % STORE IN GRID
        % --------------------------
        cd0_sum_grid(j,k)   = cd0_sum;
        cd0_misc_grid(j,k)  = CD0_misc;
        cd0_wave_grid(j,k)  = CD0_wave;
        cd0_lp_grid(j,k)    = cd0_lp;
        cd0_total_grid(j,k) = cd0_total;

    end


%% -------------------------
% STORE RESULTS IN AIRCRAFT
% --------------------------
Aircraft.constants.alt_grid_ft    = alt_ft;
Aircraft.constants.mach_grid      = mach;

Aircraft.constants.cd0_sum_grid   = cd0_sum_grid;
Aircraft.constants.cd0_misc_grid  = cd0_misc_grid;
Aircraft.constants.cd0_wave_grid  = cd0_wave_grid;
Aircraft.constants.cd0_lp_grid    = cd0_lp_grid;
Aircraft.constants.cd0_total_grid = cd0_total_grid;

%% -------------------------
% BUILD INTERPOLANT
% --------------------------
Aircraft.aero.cd0_interp = griddedInterpolant( ...
    {Aircraft.constants.alt_grid_ft, Aircraft.constants.mach_grid}, ...
    Aircraft.constants.cd0_total_grid, ...
    'linear', ...
    'nearest');

% Callable function:
% cd0_val = Aircraft.aero.cd0(alt_ft, mach)
Aircraft.aero.cd0 = @(alt_ft_query, mach_query) ...
    Aircraft.aero.cd0_interp(alt_ft_query, mach_query);

%% -------------------------
% OPTIONAL PLOTS
% --------------------------
figure;
hold on;
grid on;

plot_altitudes_ft = [0 10000 20000 30000 40000 50000];

for n = 1:length(plot_altitudes_ft)
    idx_alt = find(alt_ft == plot_altitudes_ft(n), 1);
    if ~isempty(idx_alt)
        plot(mach, cd0_total_grid(idx_alt,:), 'LineWidth', 1.6);
    end
end

xlabel('Mach number, M');
ylabel('Zero-lift drag coefficient, C_{D0}');
title('C_{D0} vs Mach at Selected Altitudes');
legend("0 ft","10,000 ft","20,000 ft","30,000 ft","40,000 ft","50,000 ft", ...
       'Location','best');

figure;
surf(mach, alt_ft, cd0_total_grid, 'EdgeColor', 'none');
xlabel('Mach number, M');
ylabel('Altitude, h [ft]');
zlabel('C_{D0}');
title('C_{D0}(h,M)');
colorbar;
view(135,30);

end

%% =========================
% LOCAL FUNCTION: SUTHERLAND VISCOSITY
% =========================
function mu = sutherland_mu(T)
% Computes dynamic viscosity of air using Sutherland's law
%
% Input:
%   T   = temperature [K]
%
% Output:
%   mu  = dynamic viscosity [kg/(m*s)]

mu_ref = 1.716e-5;   % [kg/(m*s)]
T_ref  = 273.15;     % [K]
S      = 110.4;      % [K]

mu = mu_ref * (T / T_ref)^(3/2) * (T_ref + S) / (T + S);
end
end 