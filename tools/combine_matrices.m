function [M] = combine_matrices(y, Mij)
M = y(1) * Mij{1};
for ii = 2:length(Mij)
    M = M + y(ii) * Mij{ii};
end
end