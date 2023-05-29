function [M] = sparse_grid_margin(C)
    % Input: multi-index set C
    % Ouput: reduced margin RM and margin M

    n_dim = size(C,2);

    % Create margin
    basis_vectors = eye(n_dim);
    M = zeros(0,n_dim);
    for ii = 1:size(C,1)
        M = [M ; C(ii,:) + basis_vectors];
    end
    M = setdiff(M,C, 'rows');
    M = unique(M,'rows'); % is this necesary?
end
