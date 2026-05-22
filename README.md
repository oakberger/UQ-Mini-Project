# UQ Mini Project — River Discharge Uncertainty Quantification

ETH Zürich | Uncertainty Quantification in Engineering | Spring 2026

## Project Overview

This project implements a Polynomial Chaos Expansion (PCE) surrogate model
for a river discharge formula. The goal is to quantify how uncertainty in
the channel geometry and flow properties propagates to uncertainty in the
discharge output Q.

## Model — Manning-Strickler Equation
```
Q = A · k_st · R_hy^(2/3) · J_E^(1/2)
```
For a trapezoidal cross-section:
```
A     = (b + z·h) · h 
P_wet = b + 2·h·√(1 + z²)
R_hy  = A / P_wet
```

Where:
- Q     : Discharge [m³/s]                         ← output
- A     : Cross-sectional flow area [m²]
- R_hy  : Hydraulic radius [m]
- P_wet : Wetted perimeter [m]

## Uncertain Input Parameters

| # | Symbol | Description                  | Unit        | Distribution      | Range / Parameters  |
|---|--------|------------------------------|-------------|-------------------|---------------------|
| 1 | k_st   | Strickler roughness coeff.   | m^(1/3)/s   | Normal N(μ, σ²)   | μ=35, σ=5           |
| 2 | b      | Channel bottom width         | m           | Uniform U([a, b]) | [8, 14]             |
| 3 | h      | Water depth                  | m           | Uniform U([a, b]) | [0.5, 2.5]          |
| 4 | z      | Side slope (horizontal:vert) | —           | Uniform U([a, b]) | [1.0, 2.5]          |
| 5 | J_E    | Energy slope                 | —           | Uniform U([a, b]) | [0.0005, 0.003]     |

### Parameter Justification

- **k_st**: Natural channel roughness is inherently variable and best described
  by a Normal distribution centered around a typical earth-channel value.
- **b, h, z**: Geometric properties are assumed equally likely within plausible
  ranges for a medium-sized Swiss Mittelland river — Uniform distributions apply.
- **J_E**: Channel slope varies within a narrow range typical of gentle valley
  rivers in the Swiss Mittelland — Uniform distribution applies.

## Project Structure
```text
UQ_Mini_Project/
│
├── model/
│   └── riverDischarge.m       % Manning-Strickler model
│
├── sampling/
│   ├── sampleED.m             % experimental design sampling
│   └── isoTransform.m         % isoprobabilistic transform
│
├── pce/
│   ├── buildPsi.m             % regression matrix
│   ├── computeCoeffs.m        % OLS coefficient solver
│   └── evalPCE.m              % PCE evaluation on new points
│
├── validation/
│   └── computeErrors.m        % empirical, LOO, validation error
│
├── helpers/                   % provided helper functions
│   ├── create_alphas.m
│   ├── eval_legendre.m
│   └── eval_hermite.m
│
└── experiments/
    ├── main.mlx               % master script
    ├── exp_YY_plots.mlx       % Y-Y plots
    ├── exp_errors_vs_p.mlx    % errors vs polynomial degree
    └── exp_errors_vs_n.mlx    % errors vs sample size
```
## PCE Method Summary

1. Sample experimental design using Latin Hypercube Sampling (LHS)
2. Transform inputs to standard domain via isoprobabilistic transform
3. Build multivariate polynomial basis using total degree truncation
4. Compute PCE coefficients via Ordinary Least Squares (OLS)
5. Validate surrogate using empirical, LOO, and validation errors
6. Extract mean and variance analytically from PCE coefficients

## Requirements

- MATLAB (R2021a or later recommended)
- Statistics and Machine Learning Toolbox (for lhsdesign, norminv)

## Authors

Group 10 | ETH Zürich | 2026