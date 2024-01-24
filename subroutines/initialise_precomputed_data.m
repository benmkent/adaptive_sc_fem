function precompute = initialise_precomputed_data(I_star, I, params)
%INITIALISE_PRECOMPUTED_DATA Precomputes data relating to Lagrange
%polynomials
% Inputs
%   I_star  Enhanced sparse grid
%   I       Sparse grid
%   params  Approximation parameter structure
% Outputs
%   precompute  Precomputed data

[precompute.lagrange_product_integrals, precompute.lp1] = precompute_lagrange_integrals(params); %params.test_level, params.knot_fn, params.lev2knots);
precompute = compute_lagrange_norms(I_star, I, precompute);

