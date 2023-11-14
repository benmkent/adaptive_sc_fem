function [precompute] = compute_lagrange_norms(I_star, I, precompute)
%COMPUTE_LAGRANGE_NORMS Compute the $L^2_\rho(\Gamma)$ norms of the
%multi-variate Lagrange interpolation polynomials
%
% Inputs
%   I_star      Enhanced grid points
%   I           Grid points
%   precompute  Structure of the precomputate 1D norms
% Outputs
%   precompute  Structure of the precomputed norms

lpi = precompute.lagrange_product_integrals;
I_star_r = reduce_sparse_grid(I_star);
[I_star_r] = sparse_grid_map_to_one_d_polynomials(I_star,I_star_r);
I_r = reduce_sparse_grid(I);
[I_r] = sparse_grid_map_to_one_d_polynomials(I,I_r);

%% Compute L^2 norms
% We use expansion
% (sum u_i L_i)^2 = sum(u_i^2 L_i^2) + 2 * sum_{i < j} u_i u_j L_i L_j
% and compute efficiently using our precomputed values.
n_dim = size(I_r.knots,1);
fprintf('Compute the L2 norms for I\n');
for ii = 1:I_r.size
    uv = zeros(I_r.size);
    uv(ii,ii) = 1;
    l2_norms_I(ii,1) = calc_l2rho_quick_precalc(I_r, uv, lpi, ii);
end

fprintf('Compute the L2 norm for I_star\n');
for ii = 1:I_star_r.size
    uv = zeros(I_star_r.size);
    uv(ii,ii) = 1;
    l2_norms_I_star(ii,1) = calc_l2rho_quick_precalc(I_star_r, uv, lpi, ii);
    fprintf('Lagrange polynomial L2 norm %d of %d\n',ii, I_star_r.size);
end
if I_star_r.size == I_star.size
    pts_in_both_grids_I_star = 1:I_star_r.size;
    pts_in_both_grids_I= 1:I_star_r.size;
else
    [~,pts_in_both_grids_I_star,pts_in_both_grids_I,~] = compare_sparse_grids(I_star,I_star_r,I,I_r);
end

for ii = 1:I_r.size
    u1Qv1 = zeros(I_star_r.size);
    u1Qv1(pts_in_both_grids_I_star(ii),pts_in_both_grids_I_star(ii)) = 1;
    u2Qv2 = zeros(I_r.size);
    u2Qv2(pts_in_both_grids_I(ii),pts_in_both_grids_I(ii)) = 1;
    u1Qv2 = zeros(I_star_r.size, I_r.size);
    u1Qv2(pts_in_both_grids_I_star(ii),pts_in_both_grids_I(ii)) = 1;

    interp_poly_idx = [pts_in_both_grids_I_star(ii),pts_in_both_grids_I(ii)];
    l2_norms_diff(ii,1) = calc_l2rho_diff_quick_precalc(I_star_r, I_r, u1Qv1, u2Qv2, u1Qv2, lpi, interp_poly_idx);
end

%% Approximately compute L1 norms
% This could be improved, it is not exact.
fprintf('Compute the L1 norms for I\n');
for ii = 1:I_r.size
    l1_norms_I(ii) = abs(I_r.weights(ii));
end

%% Store in precompute structure
precompute.l2_norms_I = l2_norms_I;
precompute.l2_norms_I_star = l2_norms_I_star;
precompute.l2_norms_diff = l2_norms_diff;
precompute.l1_norms_I = l1_norms_I;
end


