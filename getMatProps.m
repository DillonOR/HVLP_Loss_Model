function [perm, B_peak, rho, curie_temp, Kc, alpha, beta] = getMatProps(targetModel)
    % 1. Load data
    filename = 'magmat_database.xlsx'; 
    T = readtable(filename, 'ReadVariableNames', true); 

    % 2. Find the row based on the input 'targetModel'
    materialrow = contains(T{:, 2}, targetModel, 'IgnoreCase', true); % This finds the row I am interested in
    assignrow = find(materialrow, 1, 'first'); % This finds and gives the row a value. 

    if ~isempty(assignrow)
        % 3. Extract the specific columns (indices based on your Excel layout)
        perm       = T{assignrow, 4}; 
        B_peak     = T{assignrow, 5};
        rho        = T{assignrow, 6};
        curie_temp = T{assignrow, 7};
        Kc         = T{assignrow, 8};
        alpha      = T{assignrow, 9};
        beta       = T{assignrow, 10};
    else
        error('Material "%s" was not found in %s.', targetModel, filename);
    end
end
