function precompute = precompute_DuDv(Z_I_star, fem, precompute)    
%% Precompute L2 norms over space for all pairs of knots
    n_knots = length(Z_I_star);
%     uQv = zeros(n_knots,n_knots);
%     for ii = 1:n_knots
% %     u_z_tplusdt = [Z_I_star.u_z_tplusdt];
%     u_z_tplusdt(:,ii) = [Z_I_star{ii}.u_z_tplusdt];
%     end
%     % Why did I not include boundary???
% %     notbound = Z_I_star{1}.notbound;
% %     u_z_tplusdt = u_z_tplusdt(notbound,:);
%     Q = Z_I_star{1}.Q;
%     for ii = 1:n_knots
%         u_z_tplusdt_ii = u_z_tplusdt(:,ii);
%         for jj = 1:n_knots
%             uQv(ii,jj) = u_z_tplusdt_ii' * Q * u_z_tplusdt(:,jj);
%         end
%     end
%     precompute.uQv = uQv;

    nDiff = params.nDiff;
    %% For each collocation point pair
    for ii = 1:n_knots
        y_ii = [Z_I_star{ii}.z];
        y_diff_ii = [1; y_ii((nYConv+1):(nYConv+nYDiff))];
%         u_z_tplusdt_ii = u_z_tplusdt(:,ii);
        for jj = 1:n_knots
            y_jj = [1; Z_I_star{jj}.z];
            y_diff_jj = [1; y_jj((nYConv+1):(nYConv+nYDiff))];

%             u_z_tplusdt_jj = u_z_tplusdt(:,jj);

            %% Expand in diffusion matrices
            DD_iijj = 0;
            for kk=1:(nDiff+1)
                for ll=1:(nDiff+1)
                    DD_iijj = DD_iijj + y_diff_ii(kk) * y_diff_jj(ll) * fem.a_pars{kk,ll};
                end
            end
            DyDy(ii,jj) = DD_iijj;
%             DuDv(ii,jj) = u_z_tplusdt_ii' * DD_iijj  * u_z_tplusdt_ii;
        end
    end
    precompute.DyDy = DyDy;
end
