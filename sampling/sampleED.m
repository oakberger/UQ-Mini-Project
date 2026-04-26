    % Input:  n = sample size, params = distribution info
    % Output: X_ED (n x M) - sampled input matrix
    %         Y_ED (n x 1) - model responses
    % Does:   LHS sampling + transform + model evaluation

    function [X_ED, Y_ED] = sampleED(n, params)

M = length(params);

% 1) LHS sampling in unit hypercube
U = lhsdesign(n, M);

% 2) Transform to physical space
X_ED = zeros(n, M);

for i = 1:M
% Gibt fuer die uniform verteilten Variabeln einen Wert zwischen a und b aus
    if strcmp(params(i).type, 'uniform')
        a = params(i).a;
        b = params(i).b;
        X_ED(:,i) = a + (b - a) * U(:,i);
% Gibt fuer die normalverteilten Variabeln einen Wert gemäss Verteilung aus        
    elseif strcmp(params(i).type, 'normal')
        mu = params(i).mu;
        sigma = params(i).sigma;
        X_ED(:,i) = norminv(U(:,i), mu, sigma);
        
    else
        error('Unknown distribution type');
    end
    
end

% 3) Evaluate model
Y_ED = riverDischarge(X_ED);

end
