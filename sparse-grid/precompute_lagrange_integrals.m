function [lpi, lp1] = precompute_lagrange_integrals(params)
%PRECOMPUTE_LAGRANGE INTEGRALS Computes norms of 1D Lagrange polynomial
%products.
%
% Inputs    params  Strucutre of approximation parameters
% Outputs   lpi     Lagrange Polynomial product integral norms
%           lp1     Lagrange polynomials in one dimension L1 norms 

% Test level is the maxium level to precompute for
test_level = params.test_level;

% Generate one dimensional sets of points and weights at different levels.
knots = cell(test_level,1);
weights=cell(test_level,1);
for ii = 1:test_level
    [knots{ii}, weights{ii}] = params.knot_fn(params.lev2knots(ii)); %,-1,1,'prob');
end

% Compute grid points to to integrate with
[knots_eval, weights_eval] = params.knot_fn(test_level+1);%,-1,1,'prob');

% Evaluate each 1D polynomial on knots_eval.
Y = cell(test_level,1);
Y{1} = ones(length(knots_eval),1);
for ii=2:test_level
    x_i = knots{ii};
    Y_ii = zeros(length(knots_eval),length(x_i));
    for jj = 1:length(x_i)
        f = zeros(1,length(x_i));
        f(jj) = 1;
        % Evaluate $L_jj^{Z^ii}(y)$ at knots_eval
        Y_ii(:,jj) = modified_lagrange_formula(knots{ii}, f, knots_eval);
    end
    Y{ii} = Y_ii;
end

%% Compute 1D integrals of absolute value
% Quadrature for L^1(\Gamma) norms of each 1D polynomial.
for ii = 1:test_level
    ptwise_polys = Y{ii};
    for jj = 1:size(ptwise_polys,2);
        lp1(ii,jj) = sum(abs(ptwise_polys(:,jj)) .* weights_eval(:));
    end
end

%% Compute integrals as pairs (L2)
% Quadrature to compute $\int_{\Gamma}L_mm^{Z^ii}(y) L_{nn}^{X^jj}(y)
% \rho(y) dy$ for each possible pair of polynomials.
Z = nan(test_level,test_level,length(knots_eval),length(knots_eval));
for ii = 1:test_level
    for jj=1:test_level
        Y_ii = Y{ii};
        Y_jj = Y{jj};
        Z_ij =zeros(length(Y_ii),length(Y_jj));
        for mm = 1:size(Y_ii,2)
            for nn = 1:size(Y_jj,2);
                Z_ij(mm,nn) =(weights_eval) * (Y_ii(:,mm) .* Y_jj(:,nn));
            end
        end

        Z(ii,jj,1:size(Z_ij,1),1:size(Z_ij,2)) = Z_ij;
    end
end
lpi = Z;
end