function [S1minusS2, S1minusS2_r] = sparse_grid_subtract(S1, S2);
Sr1 = reduce_sparse_grid(S1);
Sr2 = reduce_sparse_grid(S2);

%% First check S2 is subgrid of S1
C1 = get_mi_set(S1);
C2 = get_mi_set(S2);
if isempty(setdiff(C1,C2,'rows'));
    pts_in_S2_only = [];
else
    [pts_in_S1_only,pts_in_both_grids_S1,pts_in_both_grids_S2,pts_in_S2_only] = compare_sparse_grids(S1,Sr1, S2, Sr2);
    if ~isempty(pts_in_S2_only)
        error('S2 must be subgrid of S1');
    end
end

S1minusS2 = S1;

% We now subtract the coeffs of grids in S2 from S1
jj = 1;
for ii = 1:length(S2)
    id2 = S2(ii).idx;
    while any(id2 ~= S1(jj).idx) && jj < length(S1)
        jj = jj+1;
    end
    if jj == length(S1) && any(id2 ~= S1(jj).idx);
        % Grid not present so append with negative coeff
        S2_ii_neg = S2(ii);
        S2_ii_neg.coeff = - S2_ii_neg.coeff;
        S2_ii_neg.weights = - S2_ii_neg.weights;
        S1minusS2 = [S1minusS2, S2_ii_neg];
        jj = 1;
        continue;
    else
        S1minusS2(jj).coeff = S1minusS2(jj).coeff  - S2(ii).coeff;
        S1minusS2(jj).weights = S1minusS2(jj).weights  - S2(ii).weights;
        jj = jj+1;
        continue;
    end
end

% Sort the idxs
[kk_ordered,i]=mysortrows(reshape([S1minusS2.idx],length(S1minusS2(1).idx),[])',1e-10);
S1minusS2 = S1minusS2(i);

S1minusS2([S1minusS2.coeff] ==0) = [];
S1minusS2_r = reduce_sparse_grid(S1minusS2);

%% We need to know how to map the points that remain in S1minusS2 to the points in S1
if size(S1minusS2_r.knots,2) == size(Sr1.knots,2)
    pts_in_both_grids_S1 = 1:size(S1minusS2_r.knots,2);
else
    [pts_in_S1_only,pts_in_both_grids_S1,pts_in_both_grids_S1mS2,pts_in_S1mS2_only] = compare_sparse_grids(S1,Sr1, S1minusS2, S1minusS2_r);
end
mapS2inS1 = pts_in_both_grids_S1;

% [S1minusS2_r] = MapSrToOneDPolys(S1minusS2,S1minusS2_r);
S1minusS2_r.n = mapS2inS1(S1minusS2_r.n);