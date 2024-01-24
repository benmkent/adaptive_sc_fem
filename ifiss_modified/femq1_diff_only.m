function [a, h] = femq1_diff_only(xy,ev, a_fn)
%femq1_diff_only Modified formulatuion of FEMQ1_DIFF IFISS function to
%formualte diffusion matrix only using diffusion tensor to allow 
% anisotropic diffusion.
%
% FEMQ1_DIFF vectorized bilinear coefficient matrix generator
%   [A,Q,f] = femq1_diff(xy,ev);
%   input
%          xy         vertex coordinate vector
%          ev         element mapping matrix
%   output
%          A          stiffness matrix
%          Q          mass matrix
%          f          rhs vector
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
fprintf('setting up Q1 diffusion matrices...  ')
%
% initialise global matrices
a = sparse(nvtx,nvtx);
h = sparse(nvtx,nvtx);
%
% set up 2x2 Gauss points
gpt=1.0e0/sqrt(3.0e0);
s(1) = -gpt;  t(1) = -gpt;
s(2) =  gpt;  t(2) = -gpt;
s(3) =  gpt;  t(3) =  gpt;
s(4) = -gpt;  t(4) =  gpt;
%
% inner loop over elements
for ivtx = 1:4
    xl_v(:,ivtx) = x(ev(:,ivtx));
    yl_v(:,ivtx) = y(ev(:,ivtx));
end
ae = zeros(nel,4,4);
he = zeros(nel,4,4);
%  loop over 2x2 Gauss points
for igpt = 1:4
    sigpt=s(igpt);
    tigpt=t(igpt);
    %  evaluate derivatives etc
    [jac,invjac,phi,dphidx,dphidy] = deriv(sigpt,tigpt,xl_v,yl_v);

    nel=length(xl_v(:,1));
    zero_v = zeros(nel,1); xx=zero_v; yy=xx;
    [phi_e,dphids,dphidt] = shape(sigpt,tigpt);
    for ivtx=1:4
        xx = xx + phi_e(ivtx) * xl_v(:,ivtx);
        yy = yy + phi_e(ivtx) * yl_v(:,ivtx);
    end
    
    for ii = 1:nel
        diff_mat = a_fn(xx(ii),yy(ii));

        for j = 1:4
            for i = 1:4
                ae(ii,i,j) = ae(ii,i,j)  + diff_mat(1) .* dphidx(ii,i).*dphidx(ii,j) .* invjac(ii);
                ae(ii,i,j) = ae(ii,i,j)  + diff_mat(2) .* dphidx(ii,i).*dphidy(ii,j) .* invjac(ii);
                ae(ii,i,j) = ae(ii,i,j)  + diff_mat(3) .* dphidy(ii,i).*dphidx(ii,j) .* invjac(ii);
                ae(ii,i,j) = ae(ii,i,j)  + diff_mat(4) .* dphidy(ii,i).*dphidy(ii,j) .* invjac(ii);
                he(ii,i,j) = he(ii,i,j)  + dphidx(ii,i) .* dphidx(ii,j) .* invjac(ii);
                he(ii,i,j) = he(ii,i,j)  + dphidy(ii,i) .* dphidy(ii,j) .* invjac(ii);
            end
        end
    end
    % end of Gauss point loop
end
% perform assembly of global matrix  and source vector
for krow=1:4
    nrow=ev(:,krow);
    for kcol=1:4
        ncol=ev(:,kcol);
        a = a + sparse(nrow,ncol,ae(:,krow,kcol),nvtx,nvtx);
        h = h + sparse(nrow,ncol,he(:,krow,kcol),nvtx,nvtx);
    end
end
%
return

