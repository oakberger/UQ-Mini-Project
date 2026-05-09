function [err_emp, err_loo, err_val] = computeErrors(Psi, Y_ED, c, Y_val, Y_PCE_val)
% [err_emp, err_loo, err_val] = COMPUTEERRORS(Psi, Y_ED, c, Y_val, Y_PCE_val)
% Compute three error measures for the PCE.
%
% Inputs:
%   Psi       - (n x P) regression matrix on the ED
%   Y_ED      - (n x 1) true model responses on the ED
%   c         - (P x 1) PCE coefficients from computeCoeffs
%   Y_val     - (n_val x 1) true responses on the validation set
%   Y_PCE_val - (n_val x 1) PCE responses on the validation set
%
% Outputs:
%   err_emp - relative empirical error  (optimistic: uses training data)
%   err_loo - relative leave-one-out error (better predictor of true error)
%   err_val - relative validation error (most reliable, needs separate data)

% Denominator: variance of the ED responses (normalization factor)
var_Y_ED  = sum((Y_ED  - mean(Y_ED)).^2);
var_Y_val = sum((Y_val - mean(Y_val)).^2);

% --- Empirical error ---
% PCE predictions on the training data itself -> optimistic (overfits)
Y_hat_ED = Psi * c;
err_emp  = sum((Y_ED - Y_hat_ED).^2) / var_Y_ED;

% --- LOO error (analytical, no re-fitting required) ---
% Hat matrix projects Y onto the column space of Psi
H   = Psi * ((Psi' * Psi) \ Psi');
h   = diag(H);              % leverage values (diagonal of hat matrix)
res = Y_ED - Y_hat_ED;      % residuals on ED

% Corrected residuals: removing one point changes prediction by res/(1-h)
err_loo = sum((res ./ (1 - h)).^2) / var_Y_ED;

% --- Validation error ---
% Evaluated on a completely separate set: most honest error estimate
err_val = sum((Y_val - Y_PCE_val).^2) / var_Y_val;
end