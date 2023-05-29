function [IdxSet] = GetMISet(S)
    nDims = size(S(1).knots,1);
    IdxSet = reshape([S.idx],[nDims,length(S)])';
    [~,IdxSet] = check_set_admissibility(IdxSet);
end