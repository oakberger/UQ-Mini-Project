function Q = riverDischarge(X)
    % Input:  n×5 matrix [k_st, b, h, z, J_E]
    k_st = X(:,1);
    b    = X(:,2);
    h    = X(:,3);
    z    = X(:,4);
    J_E  = X(:,5);
    
    % Output: n×1 vector of discharge values

    % Geometrie
A = h .* (b + z .* h);
U = b + 2 .* h .* sqrt(1 + z.^2);
R = A ./ U;

% Manning-Strickler
Q = A .* k_st .* (R.^(2/3)) .* sqrt(J_E);

end







