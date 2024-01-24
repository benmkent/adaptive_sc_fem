function w_eval = evaluate_w_lincomb(x1,x2,y,w_fns)
%EVALUATE_W_LINCOMB Evaluate the affine advection field
% Inputs
%   x1      vector of coordinates in x1
%   x2      vector of coordinates in x2
%   y       parameter vector
%   w_fns   cell array of advection fields
% Outputs
%   w_eval  wind field evalauted at (x1,x2) with parameter vector y

%% Compute each a_fn at the points x1,x2 for a given y
for ii = 1 : length(w_fns)
    w_fn_ii = w_fns{ii};
    w_eval_ii(ii,:,:) = y(ii)*w_fn_ii(x1,x2);
end
w_eval = squeeze(sum(w_eval_ii,1));
end