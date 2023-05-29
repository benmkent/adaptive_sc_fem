function precompute = precompute_uQv(Z_I_star, precompute)    
%% Precompute L2 norms over space for all pairs of knots
    n_knots = length(Z_I_star);
    uQv = zeros(n_knots,n_knots);
    for ii = 1:n_knots
%     u_z_tplusdt = [Z_I_star.u_z_tplusdt];
    u_z_tplusdt(:,ii) = [Z_I_star{ii}.u_z_tplusdt];
    end
    % Why did I not include boundary???
%     notbound = Z_I_star{1}.notbound;
%     u_z_tplusdt = u_z_tplusdt(notbound,:);
    Q = Z_I_star{1}.Q;
    for ii = 1:n_knots
        u_z_tplusdt_ii = u_z_tplusdt(:,ii);
        for jj = 1:n_knots
            uQv(ii,jj) = u_z_tplusdt_ii' * Q * u_z_tplusdt(:,jj);
        end
    end
    precompute.uQv = uQv;
end
