function Psi = buildPsi(Xi, alphas, params)
% Psi = BUILDPSI(Xi, alphas, params)
% Build the (n x P) regression matrix.
%
% Inputs:
%   Xi     - (n x M) standardized samples from isoTransform
%   alphas - (P x M) multi-index matrix from create_alphas(M, p)
%   params - distribution struct with fields: M, dist, a, b, mu, sigma
%
% Output:
%   Psi - (n x P) regression matrix where Psi(i,j) = Psi_j(xi^(i))

n = size(Xi, 1);
P = size(alphas, 1);
Psi = ones(n, P);   % initialize to 1 for multiplicative tensor product

for j = 1:P             % loop over basis elements (columns)
    for k = 1:params.M  % loop over input dimensions
        degree = alphas(j, k);

        if strcmp(params(k).type, 'uniform')
            % Legendre: orthonormal w.r.t. U([-1,1])
            poly_vals = eval_legendre(Xi(:, k), degree);

        elseif strcmp(params(k).type, 'normal')
            % Hermite: orthonormal w.r.t. N(0,1)
            poly_vals = eval_hermite(Xi(:, k), degree);

        else
            error('Unknown distribution type for variable %d', k);
        end

        % Tensor product across dimensions: Psi_alpha = prod_k pi_alpha_k
        Psi(:, j) = Psi(:, j) .* poly_vals;
    end
end
end