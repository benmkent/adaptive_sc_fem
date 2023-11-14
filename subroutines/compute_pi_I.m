function [pi_I] = compute_pi_I(I_star, I, precompute)
%COMPUTE_PI_I Computes interpolation error estimator
%
%   Inputs  I_star      Enhanced sparse grid
%           I           Sparse grid
%           precompute  Precomputed data

%% Preprocess enhanced sparse grid
I_star_r = reduce_sparse_grid(I_star);
[I_star_r] = sparse_grid_map_to_one_d_polynomials(I_star,I_star_r);

%% Compute interpolation error estimator
% This utilises the precomputed function u^* - u in the precompute
% dUdQdU_full precomputed data.
pi_I = calc_l2rho_quick_precalc(I_star_r,precompute.dUdQdU_full, precompute.lagrange_product_integrals);

%% Catch imaginary interpolation error estimate
% this can occur if error is almost zero as computed square norm is not
% guarenteed to be non-negative due to numerical error.
if imag(pi_I)~=0
    if imag(pi_I) < 1e-7 % Is essentially sqrt(0 - eps)
        pi_I = 0;
    else
        warning('Complex interpolation error estimate, magnitude %f',abs(pi_I));
    end
end
end
