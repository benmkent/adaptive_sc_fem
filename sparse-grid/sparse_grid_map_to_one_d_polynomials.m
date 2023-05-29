function [Sr] = sparse_grid_map_to_one_d_polynomials(S,Sr)

%% We use expansion
% (sum u_i L_i)^2 = sum(u_i^2 L_i^2) + 2 * sum_{i < j} u_i u_j L_i L_j

% First collect all Lagrange polys with IDs and values
[n_dim, n_knots] = size(Sr.knots);

inds_level = zeros(n_knots,n_dim);
inds_ind = zeros(n_knots,n_dim);
weights_all = zeros(n_knots,1);
knot_counter = 0;

for ii = 1:length(S)
    S_ii = S(ii);
    S_ii_idx = S_ii.idx;
    S_ii_m = S_ii.m;
    n_knots_ii = size(S(ii).knots,2);
    % For each idx get the corresponding list of interpolation polynomials
    % with knots
    inds = 1:S_ii_m(1);
    inds = inds(:);
    for jj = 2:n_dim
        n = size(inds,1);
        inds_old = repmat(inds,[S_ii_m(jj),1]);
        inds_new = repelem(1:S_ii_m(jj),n)';
        inds = [inds_old, inds_new];
    end
    weights_all((knot_counter+1):(knot_counter+n_knots_ii)) = S(ii).coeff;
    inds_level((knot_counter+1):(knot_counter+n_knots_ii),:) = repmat(S_ii_idx,[size(inds,1),1]);
    inds_ind((knot_counter+1):(knot_counter+n_knots_ii),:) = inds;
    knot_counter = knot_counter+n_knots_ii;
end

Sr.weights_all = weights_all;
% Sr.L_prod = L_prod;
% Sr.knots_map = knots_map;
Sr.inds_level = inds_level;
Sr.inds_ind = inds_ind;
        