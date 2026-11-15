function [Core,Ac, lc, Wa, Ap, Vc, MLT, Al, G] = Geometry(Ap_req_cm4, sheetName)
% Geometry:
%   Mode A: Geometry(Ap_req_cm4, sheetName) where Ap_req_cm4 is numeric
%           -> scans table and selects a core with Ap >= Ap_req_cm4
%
%   Mode B: Geometry(coreKey, sheetName) where coreKey is char/string
%           -> forces exact match to that core
%
% Outputs:
%   Core : selected core key/name (string)
%   Ac   : m^2
%   lc   : cm   (keep consistent with your existing downstream math)
%   Wa   : cm^2
%   Ap   : cm^4
%   Vc   : m^3
%   MLT  : cm
%   Al   : as stored in sheet (often nH/N^2)

    filename = 'magmat_database.xlsx';
    T = readtable(filename, 'Sheet', sheetName, 'ReadVariableNames', true);

    coreNames = string(T{:,1});

    if isnumeric(Ap_req_cm4)
        % -------- Mode A: scan/select by Ap requirement --------
        Ap_col = double(T{:,9}); % Ap in cm^4 (your column 9)

        idx = (Ap_col >= Ap_req_cm4) & ~isnan(Ap_col);
        if ~any(idx)
            error('No cores satisfy Ap >= %.3g cm^4 in sheet "%s".', Ap_req_cm4, sheetName);
        end

        Tok = T(idx,:);
        Ap_ok = Ap_col(idx);

        % Ranking: smallest oversize Ap (closest above requirement)
        [~,iBest] = min(Ap_ok - Ap_req_cm4);

        row = Tok(iBest,:);
        Core = string(row{1,1});

    else
        % -------- Mode B: force specific core --------
        coreKey = string(Ap_req_cm4);

        % Prefer exact match; fall back to contains if exact fails (optional)
        isExact = strcmpi(coreNames, coreKey);
        idx = find(isExact, 1, 'first');

        if isempty(idx)
            % Optional fallback: contains match (comment out if you want strict only)
            isMatch = contains(coreNames, coreKey, 'IgnoreCase', true);
            idx = find(isMatch, 1, 'first');
        end

        if isempty(idx)
            error('Core "%s" was not found in sheet "%s".', coreKey, sheetName);
        end

        row  = T(idx,:);
        Core = string(row{1,1});
    end

    % -------- Pull geometry from the chosen row --------
    % Column mapping based on your current sheet:
    % col 4: MLT (cm)
    % col 5: lc  (cm)
    % col 7: Ac  (cm^2)
    % col 8: Wa  (cm^2)
    % col 9: Ap  (cm^4)
    % col 12: Al

    Ac_cm2 = row{1,7};
    lc     = row{1,5};   % cm
    Wa     = row{1,8};   % cm^2
    Ap     = row{1,9};   % cm^4
    MLT    = row{1,4};   % cm
    Al     = row{1,12};
    G      = row{1,13};
    % Conversions (keep outputs matching your preference)
    Ac = Ac_cm2 * 1e-4;       % m^2
    Vc = Ac * (lc * 1e-2);    % m^3 (Ac in m^2, lc in m)
end

%     filename  = 'magmat_database.xlsx';
% 
%     T = readtable(filename, 'Sheet', sheetName, 'ReadVariableNames', true);
% 
%     isMatch = contains(string(T{:,1}), string(coreKey), 'IgnoreCase', true);
%     idx = find(isMatch, 1, 'first');
% 
%     if isempty(idx)
%         error('Core "%s" was not found in sheet "%s".', coreKey, sheetName);
%     end
% 
%     Ac = (T{idx, 7})*1e-4; % Ac is given in cm^2
%     lc = T{idx, 5}; %lc is also known as MPL - Magnetic Path Length - this might become different with gap size (cm)
%     Wa = T{idx, 8}; % given in cm^2
%     Ap = T{idx, 9}; % given in cm^4
%     Vc = (Ac*lc)*1E-2; % cm^2 *cm) - I multiply the answer by 1E-6 to get in m^3
%     MLT = T{idx, 4}; % cm
%     Al = T{idx, 12}; 
    

