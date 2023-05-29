function uDDv = compute_uDDv(U_I_star, precompute)    
%% Precompute L2 norms over space for all pairs of knots
    n_knots = size(U_I_star,2);
    DyDy = precompute.DyDy;
    %% For each collocation point pair
    for ii = 1:n_knots
        u_z_tplusdt_ii = U_I_star(:,ii);
        for jj = 1:n_knots
            u_z_tplusdt_jj = U_I_star(:,jj);
%             uDDv(ii,jj) = u_z_tplusdt_ii' * DyDy{ii,jj}  * u_z_tplusdt_jj;
            
        end
    end
end
