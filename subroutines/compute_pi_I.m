function [pi_I] = compute_pi_I(I_star, I, precompute)
%     uQv = precompute.uQv;
% 
%     [I_star_minus_I, I_star_minus_I_r] = sparse_grid_subtract(I_star, I);
%     [I_star_minus_I_r] = sparse_grid_map_to_one_d_polynomials(I_star_minus_I,I_star_minus_I_r);

    I_star_r = reduce_sparse_grid(I_star);
    [I_star_r] = sparse_grid_map_to_one_d_polynomials(I_star,I_star_r);

%     pi_I_2 = calc_l2rho_quick_precalc(I_star_minus_I_r,uQv, precompute.lagrange_product_integrals);
    pi_I = calc_l2rho_quick_precalc(I_star_r,precompute.dUdQdU_full, precompute.lagrange_product_integrals);

    if imag(pi_I)~=0
        if imag(pi_I) < 1e-6 % Is essentially sqrt(0 - eps)
            pi_I = 0;
        else
            warning('Complex interpolation error estimate, magnitude %f',abs(pi_I));
        end
    end
end
