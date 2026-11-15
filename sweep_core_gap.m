function [opt_core, iCore, iGap] = sweep_core_gap(AgLg, E_pulse, B_max, Lg_range, Lm_e, N)
%
    mu_0 = 4*pi*1e-7;  % H/m

    % --- Read database ---
    filename = 'magmat_database.xlsx'; 
    T = readtable(filename, 'Sheet', 'Core_Size', 'ReadVariableNames', true);


    coreName = string(T{:,1});   % core identifier/name in column A
    Ac   = (T{:,7});               % Ac in column G (cm^2)

    % --- Sweep combinations (cm^3 matching) ---
    [Ac_grid, Lg_grid] = ndgrid((Ac), Lg_range);
    AgLg_grid = Ac_grid .* Lg_grid;     % cm^3
    B_min = B_max*0.75;
   
    disp([min(Lg_grid) max(Lg_grid)])
    Np_e            =   sqrt((Lm_e.*(Lg_grid*1E-2))./(mu_0*(Ac_grid*1E-4)));        % Calculation of the number of primary turns required
    Ns_e = Np_e.*N;
    % --- Compute B correctly (convert cm^3 -> m^3 in equation) ---
   
    B_grid_T = sqrt( (2*mu_0*E_pulse) ./ (AgLg_grid*1e-6) );

    % --- Feasibility filter ---
    ok = (B_grid_T >= B_min) & (B_grid_T <= B_max) & isfinite(B_grid_T);

    % --- Matching error (cm^3) ---
    err_cm3 = abs(AgLg_grid - AgLg);

    % --- Flatten feasible results ---
    [iCore, iGap] = find(ok);

    opt_core          = table();
    opt_core.Core     = coreName(iCore);
    opt_core.Ac       = Ac(iCore);
    opt_core.Lg_mm    = (Lg_range(iGap).').*10;
    opt_core.AgLg     = AgLg_grid(ok);
    opt_core.Np_e       = Np_e(ok);
    opt_core.Ns_e       = Ns_e(ok);
    opt_core.B_T      = B_grid_T(ok);
    opt_core.err_cm3  = err_cm3(ok);

    % Sort best match first
    opt_core = sortrows(opt_core, "err_cm3", "ascend");
end
