function [I_star, I] = initialise_sparse_grid(problem, params)
    % Create grid I
    if params.l_initial == 1
        C = ones(1,problem.n);
    else
        C=multiidx_gen(problem.n,params.rule,problem.n + params.l_initial-1,1);
    end
    [adm,C] = check_set_admissibility(C);
    C = unique(C,'rows');
    I = smolyak_grid_multiidx_set(C,params.knot_fn,params.lev2knots);
    
    RM = sparse_grid_reduced_margin(C);

    [~,C_star] = check_set_admissibility([C;RM]);
    C_star = unique(C_star,'rows');
    I_star = smolyak_grid_multiidx_set(C_star,params.knot_fn,params.lev2knots);
    
    I_r = reduce_sparse_grid(I);
    I_star_r = reduce_sparse_grid(I_star);
    fprintf('Number of pts in I %d (I_star %d )\n', I_r.size, I_star_r.size )
end