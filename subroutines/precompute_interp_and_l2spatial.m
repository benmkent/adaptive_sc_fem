function precompute = precompute_interp_and_l2spatial(Z_I_star, I_star, I, fem, precompute)
%precompute_interp_and_l2spatial Precompute the vector-matrix-vector for pairs of
%collocation points with respect to mass matrix.
%
%   Inputs      Z_I_star    Enhanced sparse grid data
%               I_star      Enhanced sparse grid
%               I           Sparse grid
%               fem         Finite element data
%               precompute  Precomputed data
%   Outputs     precompute  Updated precomputed data including pairwise
%                           mass matrix weighted products.

% First we must take the current approximation and interpolate onto the new points
I_r = reduce_sparse_grid(I);
I_star_r = reduce_sparse_grid(I_star);

% Interpolate the current sln to the new pts
if I_star_r.size == I_r.size
    pts_in_I_star_only = [];
    pts_in_both_new = 1:I_star_r.size;
else
    [pts_in_I_star_only, pts_in_both_new] = compare_sparse_grids(I_star, I_star_r, I, I_r);
end
for ii = 1:length(pts_in_both_new)
    U_I(:,ii) = [Z_I_star{pts_in_both_new(ii)}.u_z_tplusdt];
end
U_I_on_I_star_only = interpolate_on_sparse_grid(I,I_r, U_I, I_star_r.knots(:,pts_in_I_star_only));

% Now must compute L2 spatial norm for pts in I_star_only
% dUQdU(pts_in_both_new) = 0;
dU = zeros(size(U_I,1),length(pts_in_I_star_only));
dUQdU = zeros(length(pts_in_I_star_only));
for ii = 1:length(pts_in_I_star_only)
    dU(:,ii) = [Z_I_star{pts_in_I_star_only(ii)}.u_z_tplusdt] - U_I_on_I_star_only(:,ii);
    for jj = 1:ii
        dUQdU(ii,jj) = dU(:,ii).' * fem.mass * dU(:,jj);
    end
end
% use symmetry
dUQdU = dUQdU + dUQdU.' - diag(diag(dUQdU));

% Now embed this in larger matrix
dUdQdU_full = sparse(I_star_r.size,I_star_r.size);
dUdQdU_full(pts_in_I_star_only,pts_in_I_star_only) = dUQdU;

precompute.dUdQdU_full = dUdQdU_full;
precompute.new_indices = pts_in_I_star_only;