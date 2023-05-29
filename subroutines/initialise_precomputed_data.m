function precompute = initialise_precomputed_data(I_star, I, params, problem)

precompute.gamma_pts_for_integration = (2*rand(problem.n,10)-1);
precompute.space_pts_for_integration = (2*rand(10,2)-1);

[precompute.lagrange_product_integrals, precompute.lp1] = precompute_lagrange_integrals(params); %params.test_level, params.knot_fn, params.lev2knots);
precompute = compute_lagrange_norms(I_star, I, precompute);

% precompute.lagrange_product_integrals = lpi;
% precompute.lp1 = lp1;
% precompute.l2_norms_I = l2_norms_I;
% precompute.l2_norms_I_star = l2_norms_I_star;
% precompute.l2_norms_diff = l2_norms_diff;
% precompute.l1_norms_I = l1_norms_I;
