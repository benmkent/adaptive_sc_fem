function y = modified_lagrange_formula(knots, f_knots, x)
%MODIFIED_LAGRANGE_FORMULA Evaluates a Lagrange polynomial using
%barycentric formulation.
% Input     knots   knots for Lagrange polynomial
%           f_knots value at each knot
%           x       points to evaluate at
% Output    y       Evaluated interpolant.
%
% See Higham 2004 https://doi.org/10.1093/imanum/24.4.547
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