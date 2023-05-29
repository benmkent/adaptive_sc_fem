function [a_pairs] = femq1_diff_pairs(xy,ev, a_fns)
%FEMQ1_DIFF vectorized bilinear coefficient matrix generator
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

% Set up pairwise diff matrices
for ii_a1 = 1:length(a_fns)
    for ii_a2 = 1:length(a_fns)
        a_fn1 = a_fns{ii_a1};
        a_fn2 = a_fns{ii_a2};

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
                diff_mat1 = a_fn1(xx(ii),yy(ii));
                diff_mat2 = a_fn2(xx(ii),yy(ii));
        
                for j = 1:4
                    for i = 1:4
                        ae(ii,i,j) = ae(ii,i,j)  + diff_mat1(1,1) .* dphidx(ii,i) .* diff_mat2(1,1) .*  dphidx(ii,j) .* invjac(ii);
                        ae(ii,i,j) = ae(ii,i,j)  + diff_mat1(1,1) .* dphidx(ii,i) .* diff_mat2(1,2) .*  dphidy(ii,j) .* invjac(ii);
                        ae(ii,i,j) = ae(ii,i,j)  + diff_mat1(1,2) .* dphidy(ii,i) .* diff_mat2(1,1) .*  dphidx(ii,j) .* invjac(ii);
                        ae(ii,i,j) = ae(ii,i,j)  + diff_mat1(1,2) .* dphidy(ii,i) .* diff_mat2(1,2) .*  dphidy(ii,j) .* invjac(ii);
                        ae(ii,i,j) = ae(ii,i,j)  + diff_mat1(2,1) .* dphidx(ii,i) .* diff_mat2(2,1) .*  dphidx(ii,j) .* invjac(ii);
                        ae(ii,i,j) = ae(ii,i,j)  + diff_mat1(2,1) .* dphidx(ii,i) .* diff_mat2(2,2) .*  dphidy(ii,j) .* invjac(ii);
                        ae(ii,i,j) = ae(ii,i,j)  + diff_mat1(2,2) .* dphidy(ii,i) .* diff_mat2(2,1) .*  dphidx(ii,j) .* invjac(ii);
                        ae(ii,i,j) = ae(ii,i,j)  + diff_mat1(2,2) .* dphidy(ii,i) .* diff_mat2(2,2) .*  dphidy(ii,j) .* invjac(ii);
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


        %% Matrix of matrices
        a_pairs{ii_a1,ii_a2} = a; 
    end
end
%
return

