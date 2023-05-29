function J = mark_I(pi_I_alpha, work, RM, E, params)
    total_error = sum(pi_I_alpha);
    
%     cost = params.marking_fn(pi_I_alpha,work);
    cost = pi_I_alpha;
    [~,inds] = sort(cost,'descend');
    
    cumulative_error = 0;
    jj = 0;
    while (total_error - cumulative_error) > total_error * params.marking_factor
        jj = jj + 1;
        cumulative_error = sum(pi_I_alpha(inds(1:jj)));
    end
    marked = inds(1:jj);

    J = RM(marked,:);