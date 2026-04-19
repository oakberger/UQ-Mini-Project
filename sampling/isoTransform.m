function Xi = isoTransform(X, params)
    % Input:  X (n×M) in original space, params = distribution info
    % Output: Xi (n×M) in standard space [-1,1] or N(0,1)
    % Handles both uniform and normal distributions

    [n, M] = size(X);
    Xi = zeros(n, M);

    for i = 1:M
    
    if strcmp(params(i).type, 'uniform')
        a = params(i).a;
        b = params(i).b;
        Xi(:,i) = (2*X(:,i) - (a + b)) / (b - a);
        
    elseif strcmp(params(i).type, 'normal')
        mu = params(i).mu;
        sigma = params(i).sigma;
        Xi(:,i) = (X(:,i) - mu) / sigma;
        
    else
        error('Unknown distribution type');
    end
    
end
end






