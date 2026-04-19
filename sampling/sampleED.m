    % Input:  n = sample size, params = distribution info
    % Output: X_ED (n×M), Y_ED (n×1)
    % Does:   LHS sampling + transform + model evaluation

    function [X_ED, Y_ED] = sampleED(n, params)
% SAMPLEED Creates an experimental design using LHS
%
% Input:
%   n        - number of samples
%   params   - struct array with fields:
%              .type = 'uniform' or 'normal'
%              .a, .b (for uniform)
%              .mu, .sigma (for normal)
%   modelFun - function handle (e.g. @myModel)
%
% Output:
%   X_ED (n x M) - sampled input matrix
%   Y_ED (n x 1) - model responses

M = length(params);

% 1) LHS sampling in unit hypercube
U = lhsdesign(n, M);

% 2) Transform to physical space
X_ED = zeros(n, M);

for i = 1:M
    
    if strcmp(params(i).type, 'uniform')
        a = params(i).a;
        b = params(i).b;
        X_ED(:,i) = a + (b - a) * U(:,i);
        
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
