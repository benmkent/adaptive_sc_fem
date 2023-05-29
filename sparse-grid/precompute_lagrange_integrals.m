function [lpi, lp1] = precompute_lagrange_integrals(params)

test_level = params.test_level;

% Generate one dimensional sets of points
for ii = 1:test_level
    [knots{ii} weights{ii}] = params.knot_fn(params.lev2knots(ii)); %,-1,1,'prob');
end
% [knots_eval, weights_eval] = params.knot_fn(params.lev2knots(test_level+1));%,-1,1,'prob');
% [knots_eval, weights_eval] = params.knot_fn([2*length(knots{test_level})+1);%,-1,1,'prob');
[knots_eval, weights_eval] = params.knot_fn(31);%,-1,1,'prob');

% Eval on grid of test_level+1
Y{1} = ones(length(knots_eval),1);
Y_fn{1}{1} = @(y) ones(size(y));
for ii=2:test_level
    x_i = knots{ii};
    for jj = 1:length(x_i);
        f = zeros(1,length(x_i));
        f(jj) = 1;
        Y_ii(:,jj) = modified_lagrange_formula(knots{ii}, f, knots_eval);
        Y_ii_fn{jj} = @(y) modified_lagrange_formula(knots{ii}, f, y);
    end
    Y{ii} = Y_ii;
    Y_fn{ii} = Y_ii_fn;
end

%% Compute 1D integrals of absolute value
for ii = 1:test_level
    ptwise_polys = Y{ii};
    for jj = 1:size(ptwise_polys,2);
        lp1(ii,jj) = sum(abs(ptwise_polys(:,jj)) .* weights_eval(:));
    end
end

%% Compute integrals as pairs (L2)
Z = nan(test_level,test_level,length(knots_eval),length(knots_eval));
for ii = 1:test_level
    for jj=1:test_level
        Y_ii = Y{ii};
%         Y_ii_fn = Y_fn{ii};
        Y_jj = Y{jj};
%         Y_jj_fn = Y_fn{jj};
        Z_ij =zeros(length(Y_ii),length(Y_jj));
%         Z_ij_fn = [];
        for mm = 1:size(Y_ii,2)
            for nn = 1:size(Y_jj,2);
%                 Z_ij(mm,nn) =(rho(knots_eval) .* weights_eval) * (Y_ii(:,mm) .* Y_jj(:,nn));
                Z_ij(mm,nn) =(weights_eval) * (Y_ii(:,mm) .* Y_jj(:,nn));
%                 fn1 = Y_ii_fn{mm};
%                 fn2 = Y_jj_fn{nn};
%                 fn = @(y) fn1(y).*fn2(y) * 0.5;
%                 Z_ij_fn(mm,nn) = integral(fn,-1,1);
            end
        end
%         Z{ii,jj} = Z_ij;
%         Z{jj,ii} = Z_ij';
        Z(ii,jj,1:size(Z_ij,1),1:size(Z_ij,2)) = Z_ij;
%         Z(jj,ii,1:size(Z_ij,2),1:size(Z_ij,1)) = Z_ij';
%         Zalt{ii,jj} = Z_ij_fn;
%         Zalt{jj,ii} = Z_ij_fn';
    end
end
lpi = Z;
end

function y = modified_lagrange_formula(knots, f_knots, x)
%% See Higham2004
n = length(knots);
l = @(z) prod(z - knots);
w_j = @(knots,j) 1/prod(knots(j) - knots((1:n)~=j));

for ii = 1:n
    w(ii) = w_j(knots,ii);
end

for ii = 1:length(x)
    if any(x(ii) == knots)
        jj = find(x(ii) == knots,1);
        y(ii) = f_knots(jj);
    else
        y(ii) = l(x(ii)) * sum(f_knots.*w./(x(ii)-knots));
    end
end
end