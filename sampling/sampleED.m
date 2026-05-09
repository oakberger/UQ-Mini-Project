% Input:  n = sample size, params = distribution info
% Output: X_ED (n x M) - sampled input matrix
%         Y_ED (n x 1) - model responses
% Does:   LHS sampling + transform + model evaluation

function [X_ED, Y_ED] = sampleED(n, params, model)

M = numel(params);

% 1) LHS sampling in unit hypercube
U = lhsdesign(n, M);

% 2) Transform to physical space
X_ED = zeros(n, M);

for i = 1:M
    % Gibt fuer die uniform verteilten Variabeln einen Wert zwischen a und b aus
    if strcmp(params(i).type, 'uniform')
        X_ED(:,i) = params(i).a + (params(i).b - params(i).a) * U(:,i);
    elseif strcmp(params(i).type, 'normal')
        X_ED(:,i) = norminv(U(:,i), params(i).mu, params(i).sigma);
    else
        error('Unknown distribution type');
    end

end

% 3) Evaluate model
Y_ED = model(X_ED);

end
