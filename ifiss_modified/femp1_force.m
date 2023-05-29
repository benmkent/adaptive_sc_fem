function [f] = femp1_force(xy,ev, f_fn)
%FEMQ1_DIFF vectorized bilinear coefficient matrix generator
%   [A,Q,f] = femq1_diff(xy,ev);
%   input
%          xy         vertex coordinate vector
%          ev         element mapping matrix
%   output
%          f          f vector
%
%   Natural boundary conditions apply. Dirichlet conditions
%   must be explicitly enforced by calling function nonzerobc.
%   IFISS function: DJS; 4 March 2005.
% Copyright (c) 2005 D.J. Silvester, H.C. Elman, A. Ramage
x=xy(:,1); y=xy(:,2);
nvtx=length(x);
nel=length(ev(:,1));
lx=max(x)-min(x); ly=max(y)-min(y);
hx=max(diff(x)); hy=max(diff(y));
fprintf('setting up f matrices...  ')
%
% initialise global matrices
f = sparse(nvtx,1);
%
% set up Gauss points
ngpt=7; [s,t,wt]=triangular_gausspoints(ngpt);

% inner loop over elements
for ivtx = 1:3
    xl_v(:,ivtx) = x(ev(:,ivtx));
    yl_v(:,ivtx) = y(ev(:,ivtx));
end
fe = zeros(nel,3);
%  loop over 2x2 Gauss points
for igpt = 1:ngpt
    sigpt=s(igpt);
    tigpt=t(igpt);
    %  evaluate derivatives etc
    [jac,invjac,phi,dphidx,dphidy] = tderiv(sigpt,tigpt,xl_v,yl_v);

    nel=length(xl_v(:,1));
    nevtx = size(xl_v,2);
    zero_v = zeros(nel,1); xx=zero_v; yy=xx;
    [phi_e,dphids,dphidt] = tshape(sigpt,tigpt);
    for ivtx=1:nevtx
        xx = xx + phi_e(ivtx) * xl_v(:,ivtx);
        yy = yy + phi_e(ivtx) * yl_v(:,ivtx);
    end

    f_gpt = f_fn([xx,yy]);

    for i = 1:3
            fe(:,i) = fe(:,i)  + wt(igpt)*f_gpt.*phi(:,i) .* jac(:);
    end
    % end of Gauss point loop
end
% perform assembly of global matrix  and source vector
for krow=1:3
    nrow=ev(:,krow);
    f = f + sparse(nrow,1,fe(:,krow),nvtx,1);
end
%
return

