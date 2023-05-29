function [pi_I_alpha, work, RM] = compute_pi_I_alpha(I_star, I, params, precompute)
    % Get full grids
    n_dim = size(I(1).knots,1);
    I_r = reduce_sparse_grid(I);
    I_star_r = reduce_sparse_grid(I_star);
    C_star = get_mi_set(I_star);
    C = get_mi_set(I);
    [RM] = setdiff(C_star,C,'rows');
    
%     uQv = precompute.uQv;
    dUQdU = precompute.dUdQdU_full;

    for ii = 1:size(RM,1);
        RMii = RM(ii,:);
        [adm,C_ordered] = check_set_admissibility([C; RMii]);
%         if adm == false
%             error('Should be admissible');
%         end
        C_ordered = unique(C_ordered,'rows');
        I_RMii = smolyak_grid_multiidx_set(C_ordered,params.knot_fn,params.lev2knots, I);
        I_RMii_r = reduce_sparse_grid(I_RMii);
        [I_RMii_r] = sparse_grid_map_to_one_d_polynomials(I_RMii,I_RMii_r);

        if size(RM,1) ~=1
            [ ~,inds_I_RMii_in_I_star,~,~] = compare_sparse_grids(I_star,I_star_r,I_RMii,I_RMii_r);
        else
            inds_I_RMii_in_I_star = 1:I_star_r.size;
        end
        
%         u1Qv1 = uQv(inds_I_RMii_in_I_star,inds_I_RMii_in_I_star);
        dUQdU_ii = dUQdU(inds_I_RMii_in_I_star,inds_I_RMii_in_I_star);

%         [I_RMii_minus_I, I_RMii_minus_I_r] = sparse_grid_subtract(I_RMii, I);
%         [I_RMii_minus_I_r] = sparse_grid_map_to_one_d_polynomials(I_RMii_minus_I,I_RMii_minus_I_r);

%         pi_I_alpha_2(ii) = calc_l2rho_quick_precalc(I_RMii_minus_I_r,u1Qv1, precompute.lagrange_product_integrals);
        pi_I_alpha(ii) = calc_l2rho_quick_precalc(I_RMii_r,dUQdU_ii, precompute.lagrange_product_integrals);
        
        work(ii) = I_RMii_r.size - I_r.size;
    end
end
