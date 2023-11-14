function a_eval = evaluate_a_lincomb(x1,x2,y,a_fns)
%EVALUATE_A_LINCOMB Evaluates a y weighted affine combination of diffusion
%fields
%
% Inputs    x1     x1 coordinates
%           x2     x2 coordinates
%           y      parameter vector y
%           a_fns  cell array of diffusion fields
% Outputs   a_eval diffusion field at (x1,x2)

%% Compute each a_fn at the points x1,x2 for a given y
for ii = 1 : length(a_fns)
    a_fn_ii = a_fns{ii};
    %         for kk = 1 : size(x1,1)
    a_eval_ii(ii,:,:) = y(ii)* a_fn_ii(x1,x2);
    %         end
end
a_eval = squeeze(sum(a_eval_ii,1));
end