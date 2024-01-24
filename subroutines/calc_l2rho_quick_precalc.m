function [eval, eval_no_sqrt] = calc_l2rho_quick_precalc(Sr,uQv, Lagrange_Product_Integrals, varargin)
%calc_l2rho_quick_precalc Quick computation of L^2_\rho(\Gamma;\Vert \cdot
%  \Vert_{Q}^2) norms
%
% Inputs    Sr      reduced sparse grid
%           uQv     precomputed uQv vector-matrix-vector products
%           Lagrange_Product_Integrals  Precomputed integrals
%           varargin   Can input a vector of indices for collocation points
%           with non-zero solution to speed up computation.
%

% Need number of items in Sr_1, Sr_2
[n_dim, n_knots] = size(Sr.knots);
n_S = length(Sr.n);
Sr_n = Sr.n;
% Get grid point to polynomial mapping data
idx = Sr.inds_level;
ind = Sr.inds_ind;
coeffs = Sr.weights_all;

L_prod_ij = zeros(n_dim,1);

% Set up iterator over all grid points, or altertatively using the known
% non-zero indices input to the function.
iter = 1:n_S;
if nargin == 4
    iter = iter(Sr.n == varargin{1});
end

sum_L_prod = 0;

%% Iteratate over the grid points.
for ii = iter
    ind_ii = ind(ii,:);
    idx_ii = idx(ii,:);
    n_ii = Sr_n(ii);
    uQv_ii = uQv(n_ii,:);
    if all(uQv_ii == 0)
%         L_prod(ii,:) = 0;
    else
        for jj = iter
            ind_jj = ind(jj,:);
            idx_jj = idx(jj,:);
            n_jj = Sr_n(jj);
            % Select appropriate integral
    %         Z_ij = Lagrange_Product_Integrals{idx1,idx2};
            % Now select appropriate fn vales and calc weighted inner product
            % for ii we need knots Sr_1.n(ii)
            % for jj we need Sr_2.n(jj)
    %         prod_Q(ii,jj) = f_on_pts(:,Sr.n(ii))'*Q*f_on_pts(:,Sr.n(jj));
            prod_Q_ii_jj = uQv_ii(n_jj);
            % Now multiply by appropriate Lagrange poly integrals
            for kk=1:n_dim
                L_prod_ij(kk) = Lagrange_Product_Integrals(idx_ii(kk),idx_jj(kk),ind_ii(kk),ind_jj(kk));
            end
            L_prod_ii_jj = coeffs(ii) * coeffs(jj) * prod_Q_ii_jj * prod(L_prod_ij);
            sum_L_prod = sum_L_prod + L_prod_ii_jj;
        end
    end
end

%% Return the sum and the sqrt
eval_no_sqrt = sum_L_prod;
eval = sqrt(eval_no_sqrt);