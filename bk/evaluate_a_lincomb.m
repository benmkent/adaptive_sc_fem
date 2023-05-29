function a_eval = evaluate_a_lincomb(x1,x2,y,a_fns)
    %% Compute each a_fn at the points x1,x2 for a given y
    for ii = 1 : length(a_fns)
        a_fn_ii = a_fns{ii};
%         for kk = 1 : size(x1,1)
            a_eval_ii(ii,:,:) = y(ii)* a_fn_ii(x1,x2);
%         end
    end
    a_eval = squeeze(sum(a_eval_ii,1));
end