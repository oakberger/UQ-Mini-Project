function c = computeCoeffs(Psi, Y_ED)
% c = COMPUTECOEFFS(Psi, Y_ED)
% Compute PCE coefficients by Ordinary Least Squares.
%
% Inputs:
%   Psi  - (n x P) regression matrix from buildPsi
%   Y_ED - (n x 1) model responses on the experimental design
%
% Output:
%   c - (P x 1) vector of PCE coefficients
%
% Note: moments are extracted in main as:
%   mean_PCE = c(1)
%   var_PCE  = sum(c(2:end).^2)

c = Psi \ Y_ED;

end