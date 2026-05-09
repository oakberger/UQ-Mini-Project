function Xi = isoTransform(X, params)
    % Input:  X (n×M) in original space, params = distribution info
    % Output: Xi (n×M) in standard space [-1,1] or N(0,1)
    % Handles both uniform and normal distributions

    [n, M] = size(X);
    Xi = zeros(n, M);

    for i = 1:M
    
        if strcmp(params(i).type, 'uniform')
            Xi(:,i) = (2*X(:,i) - (params(i).a + params(i).b)) / (params(i).b - params(i).a);
        elseif strcmp(params(i).type, 'normal')
            Xi(:,i) = (X(:,i) - params(i).mu) / params(i).sigma;
        else
            error('Unknown distribution type');
        end
    
    end
end






