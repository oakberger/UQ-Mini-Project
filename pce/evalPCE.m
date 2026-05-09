function Y_PCE = evalPCE(X_new, c, alphas, params)
% Y_PCE = EVALPCE(X_new, c, alphas, params)
% Evaluate the trained PCE surrogate at new physical-space inputs.
%
% Inputs:
%   X_new  - (n_new x M) new samples in physical space
%   c      - (P x 1) PCE coefficients from computeCoeffs
%   alphas - (P x M) multi-index matrix from create_alphas
%   params - distribution struct (same as used everywhere else)
%
% Output:
%   Y_PCE - (n_new x 1) PCE predictions

% Step 1: transform to standard space (same transform as training)
Xi_new = isoTransform(X_new, params);

% Step 2: build regression matrix for new points
Psi_new = buildPsi(Xi_new, alphas, params);

% Step 3: evaluate surrogate
Y_PCE = Psi_new * c;
end