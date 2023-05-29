function [pi_delta, pi_I_delta, Z_I_star] = compute_pi_delta(Z_I_star, Z_I_star_lofi, I_star, I, params, precompute)
%COMPUTE_PI_DELTA Computes timestepping contribution to error
%   Detailed explanation goes here
    for ii = 1:length(Z_I_star)
        ge_error_lofi = Z_I_star{ii}.u_z_tplusdt - Z_I_star_lofi{ii}.u_z_tplusdt;
%         ge_error_lofi_notbound = ge_error_lofi(Z_I_star{ii}.notbound);
        ge_error_lofi_notbound = ge_error_lofi; % WHY DID I NOT INCLUDE BOUNDARY?
        Q = Z_I_star{ii}.Q;

        ge_error_lofi_norm(ii) = sqrt(ge_error_lofi_notbound' * Q * ge_error_lofi_notbound);
        
        p = 2;
        ge_error_norm(ii) = (params.letol/params.letol_lofi)^(p/(p+1)) * ge_error_lofi_norm(ii);
        Z_I_star{ii}.ge_estimate = ge_error_norm(ii);
    end
    
        I_r = reduce_sparse_grid(I);
        I_star_r = reduce_sparse_grid(I_star);
        [I_star_only, I_in_I_star] = compare_sparse_grids(I_star, I_star_r, I, I_r);

    if params.simplified_estimator == 0
        % Now must combine with Lagrange interpolation polynomial norms.
        pi_delta = sum(ge_error_norm(I_in_I_star) .* precompute.l2_norms_I(:)');
        pi_I_delta = sum(ge_error_norm(I_star_only) .* precompute.l2_norms_I_star(I_star_only)') + sum(ge_error_norm(I_in_I_star) .* precompute.l2_norms_diff(:)');
    else
        pi_delta = mean(ge_error_norm(I_in_I_star));
        pi_I_delta = 0;
    end

    pi_delta = abs(pi_delta);
    pi_I_delta = abs(pi_I_delta);
end