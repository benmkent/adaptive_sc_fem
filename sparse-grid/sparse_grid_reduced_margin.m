function [RM, M] = sparse_grid_reduced_margin(C)
%SPARSE_GRID_REDUCED_MARGIN Identifies reduced margin for C
% Input: multi-index set C
% Ouput: reduced margin RM and margin M

%% Create margin
n_dim = size(C,2);
basis_vectors = eye(n_dim);
M = zeros(0,n_dim);
for ii = 1:size(C,1)
    M = [M ; C(ii,:) + basis_vectors];
end
M = setdiff(M,C, 'rows');
M = unique(M,'rows'); % is this necesary?

%% Create reduced margin
for ii = 1:size(M,1)
    [adm(ii), completed_set, missing_set] = check_index_admissibility(M(ii,:),C);
end
% Admissible indices in M are RM
RM = M(adm == 1,:);
end
