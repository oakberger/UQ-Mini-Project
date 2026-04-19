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


%Kommt dieser Teil in "main"?
params(1).type = 'normal';
params(1).mu = 35;
params(1).sigma = 5;

params(2).type = 'uniform';
params(2).a = 8;
params(2).b = 14;

params(3).type = 'uniform';
params(3).a = 0.5;
params(3).b = 2.5;

params(4).type = 'uniform';
params(4).a = 1.0;
params(4).b = 2.5;

params(5).type = 'uniform';
params(5).a = 0.0005;
params(5).b = 0.003;
