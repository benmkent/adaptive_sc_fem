function [eval, prod_Q] = calc_l2rho_diff_quick_precalc(Sr_1, Sr_2,  u1Qv1, u2Qv2, u1Qv2, Lagrange_Product_Integrals, varargin)
% Optional interp_poly_idx speeds up processing by skipping zero value indices in for loops.

%% Calculate u1^2 + u2^2
if nargin == 7
    interp_poly_idx = varargin{1};
    [eval_1, eval_no_sqrt_1] = calc_l2rho_quick_precalc(Sr_1, u1Qv1, Lagrange_Product_Integrals, interp_poly_idx(1));
    [eval_2, eval_no_sqrt_2] = calc_l2rho_quick_precalc(Sr_2, u2Qv2, Lagrange_Product_Integrals, interp_poly_idx(2));
else
    [eval_1, eval_no_sqrt_1] = calc_l2rho_quick_precalc(Sr_1, u1Qv1, Lagrange_Product_Integrals);
    [eval_2, eval_no_sqrt_2] = calc_l2rho_quick_precalc(Sr_2, u2Qv2, Lagrange_Product_Integrals);
end

%% Next calculate the product entries
% Need number of items in Sr_1, Sr_2
[n_dim] = size(Sr_1.knots,1);
n_S_1 = length(Sr_1.n);
n_S_2 = length(Sr_2.n);
Sr_1_n = Sr_1.n;
Sr_2_n = Sr_2.n;

idx1 = Sr_1.inds_level;
idx2 = Sr_2.inds_level;
ind1 = Sr_1.inds_ind;
ind2 = Sr_2.inds_ind;
coeffs_1 = Sr_1.weights_all;
coeffs_2 = Sr_2.weights_all;

% prod_Q = zeros(n_S_1,n_S_2);
L_prod_ij = zeros(n_dim,1);
% L_prod = zeros(n_S_1,n_S_2);

uQv_reduced = u1Qv2;

% If general poly need to consdier all points in grid
    iter_ii = 1:n_S_1;
    iter_jj = 1:n_S_2;
if nargin == 7
    % If a interpolation poly lots of values at knots are zero so only
    % consider those that are equal to 1.
    iter_ii = iter_ii(Sr_1_n == interp_poly_idx(1));
    iter_jj = iter_jj(Sr_2_n == interp_poly_idx(2));
end

sum_L_prod = 0;

for ii = iter_ii
    for jj = iter_jj
        % Select appropriate integral
%         Z_ij = Lagrange_Product_Integrals{idx1,idx2};
        % Now select appropriate fn vales and calc weighted inner product
        % for ii we need knots Sr_1.n(ii)
        % for jj we need Sr_2.n(jj)
%         prod_Q(ii,jj) = fn_on_pts_1(:,Sr_1.n(ii))'*Q*fn_on_pts_2(:,Sr_2.n(jj));
%         prod_Q(ii,jj) = uQv_reduced(Sr_1_n(ii),Sr_2_n(jj));
        prod_Q_ii_jj = uQv_reduced(Sr_1_n(ii),Sr_2_n(jj));
        % Now multiply by appropriate Lagrange poly integrals
        for kk=1:n_dim
%             Z_ij = Lagrange_Product_Integrals{idx1(ii,kk),idx2(jj,kk)};
%             L_prod_ij(kk) = Z_ij(ind1(ii,kk),ind2(jj,kk));
            L_prod_ij(kk) = Lagrange_Product_Integrals(idx1(ii,kk),idx2(jj,kk),ind1(ii,kk),ind2(jj,kk));
        end
        L_prod_ii_jj = coeffs_1(ii) * coeffs_2(jj) * prod_Q_ii_jj * prod(L_prod_ij);
%         L_prod(ii,jj) = coeffs_1(ii) * coeffs_2(jj) * prod_Q(ii,jj) * prod(L_prod_ij);
%         sum_L_prod = sum_L_prod + L_prod(ii,jj);
        sum_L_prod = sum_L_prod + L_prod_ii_jj;
    end
end

eval = sqrt(eval_no_sqrt_1 + eval_no_sqrt_2 - 2*sum_L_prod);