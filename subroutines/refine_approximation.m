function [Z_I_star_refined, I_new ,I_star_new] = refine_approximation(Z_I_star, I, I_star, J, problem, params, fem, hifi)
    %% Define refined grid
    C = get_mi_set(I);
    C_new = [C;J];
    [~,C_new_ordered] = check_set_admissibility(C_new);
    C_new_ordered = unique(C_new_ordered,'rows');
    I_new = smolyak_grid_multiidx_set(C_new_ordered,params.knot_fn,params.lev2knots);
    
    C_rm = sparse_grid_reduced_margin(C_new_ordered);
    C_star = get_mi_set(I_star);
    [~,C_star_new_ordered] = check_set_admissibility([C_star; C_rm]);
    C_star_new_ordered = unique(C_star_new_ordered,'rows');
    I_star_new = smolyak_grid_multiidx_set(C_star_new_ordered,params.knot_fn,params.lev2knots);
    
    I_star_new_r = reduce_sparse_grid(I_star_new);
    I_star_r = reduce_sparse_grid(I_star);
    if I_star_new_r.size == I_star_r.size
        new_pts = [];
        old_pts = 1:I_star_new_r.size;
        old_pts_in_I_star = old_pts;    
    else
        [new_pts,old_pts,old_pts_in_I_star,~] = compare_sparse_grids(I_star_new, I_star_new_r, I_star, I_star_r);
    end

    %% Initialise structure for ALL points
    Z_I_star_refined = initialise_approximation_structure(problem, params, fem, I_star_new);

    % Copy current points structure into appropriate places
    for ii = 1:length(old_pts);
        Z_I_star_refined{old_pts(ii)} = Z_I_star{old_pts_in_I_star(ii)};
    end

    if params.adapt_interp == 1 && hifi == 1
        % Construct u on old grid
        t = Z_I_star{1}.t_z(1);
        for ii = 1:length(old_pts)
            U(:,ii) = Z_I_star_refined{old_pts(ii)}.u_z(:,1);
        end
        U_new = interpolate_on_sparse_grid(I_star, I_star_r,U, I_star_new_r.knots(:,new_pts));
        for ii = 1:length(new_pts)
            Z_I_star_refined{new_pts(ii)}.u_z = U_new(:,ii);
            Z_I_star_refined{new_pts(ii)}.t_z = t;
            Z_I_star_refined{new_pts(ii)}.u_z_tplusdt = U_new(:,ii);
            Z_I_star_refined{new_pts(ii)}.tplusdt = t;
        end
    end
end