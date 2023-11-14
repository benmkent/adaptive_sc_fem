function [a, h] = femp1_diff_only(xy,evt, a_fn)
%FEMP1_DIFF_ONLY set up linear anisotropic diffusion matrices for P1 mesh.
% Based upon TIFISS code.
%FEMP1_DIFF  set up linear anisotropic diffusion matrices
%   [A,Q,f,Ae,Qe] = femp1_diff(xy,evt);
%   input
%          xy         vertex coordinate vector  
%          evt        element mapping matrix
%   output
%          A          diffusion matrix
%          Q          mass matrix 
%          f          rhs vector
%          Ae         element diffusion matrices
%          Qe         element mass matrices
%
%   Natural boundary conditions apply. Dirichlet conditions
%   must be explicitly enforced by calling function nonzerobc.
%    TIFISS function: DJS; 3 March 2017.
% Copyright (c) 2007 C.E. Powell, D.J. Silvester
x=xy(:,1); y=xy(:,2);
nvtx=length(x);
nel=length(evt(:,1));
fprintf('setting up P1 diffusion matrices...  ')
%
% initialise global matrices
      a = sparse(nvtx,nvtx);
      h = sparse(nvtx,nvtx);
%
% inner loop over elements    
        for ivtx = 1:3
        xl_v(:,ivtx) = x(evt(:,ivtx));
        yl_v(:,ivtx) = y(evt(:,ivtx)); 
	    end
        ae = zeros(nel,3,3);
        he = zeros(nel,3,3);
%
%  Gauss rule integration
%          sigpt=1/3;
%          tigpt=1/3;
%          wt=1/2;
nngpt=7; [s,t,wt]=triangular_gausspoints(nngpt);

for igpt=1:nngpt
   sigpt=s(igpt);
    tigpt=t(igpt);
    wtigpt=wt(igpt);
%  evaluate derivatives etc
         [jac,invjac,phi,dphidx,dphidy] = tderiv(sigpt,tigpt,xl_v,yl_v);
%          [diffx,diffy] = tgauss_adiff(sigpt,tigpt,xl_v,yl_v); 
         nel=length(xl_v(:,1));
         zero_v = zeros(nel,1); xx=zero_v; yy=xx;
         [xi,dxids,dxidt] = tshape(sigpt,tigpt);
         for ivtx=1:3
             xx = xx + xi(ivtx) * xl_v(:,ivtx);
             yy = yy + xi(ivtx) * yl_v(:,ivtx);
         end
%         diffx = a_fn(xx,yy);
%         diffy = a_fn(xx,yy);

%         for kk = 1:size(xx,1);
            diff_matrix = a_fn(xx(:),yy(:));
             for j = 1:3
                   for i = 1:3
                   ae(:,i,j) = ae(:,i,j) + wtigpt* diff_matrix(:,1) .* dphidx(:,i).*dphidx(:,j) .* invjac(:);
                   ae(:,i,j) = ae(:,i,j) + wtigpt* diff_matrix(:,2) .* dphidx(:,i).*dphidy(:,j) .* invjac(:);
                   ae(:,i,j) = ae(:,i,j) + wtigpt* diff_matrix(:,3) .* dphidy(:,i).*dphidx(:,j) .* invjac(:);
                   ae(:,i,j) = ae(:,i,j) + wtigpt* diff_matrix(:,4) .*dphidy(:,i).*dphidy(:,j) .* invjac(:);
                   he(:,i,j) = he(:,i,j) + wtigpt* dphidx(:,i).*dphidx(:,j) .* invjac(:);
                   he(:,i,j) = he(:,i,j) + wtigpt* dphidy(:,i).*dphidy(:,j) .* invjac(:);
                   end
             end
%         end
end
%
% perform assembly of global matrix  and source vector
      for krow=1:3
	  nrow=evt(:,krow);	 
          for kcol=1:3
		  ncol=evt(:,kcol);	  
          a = a + sparse(nrow,ncol,ae(:,krow,kcol),nvtx,nvtx);
          h = h + sparse(nrow,ncol,he(:,krow,kcol),nvtx,nvtx);
          end
%      f(nrow,1) = f(nrow,1) + fe(:,krow)
      end
%
return

