function [wire_rp, wire_rs] = getwire_resistance(A_wp, A_ws)
    % 1. Load data
    filename = 'magmat_database.xlsx'; 
    T = readtable(filename, 'Sheet', 'Wire_Size', 'ReadVariableNames', true);
    
    dia = T{:,2};
    res = T{:,4};

    % 2. Find the row based on the input approx size of wire
   [~, pos_wp] = min(abs(dia-A_wp)); 
   [~, pos_ws] = min(abs(dia-A_ws)); 
    
   wire_rp = res(pos_wp);
   wire_rs = res(pos_ws);
    
end
