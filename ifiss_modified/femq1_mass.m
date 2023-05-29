function [Q] = femq1_mass(xy,ev)
%FEMQ1_DIFF vectorized bilinear coefficient matrix generator
%   [A,Q,f] = femq1_diff(xy,ev);
%   input
%          xy         vertex coordinate vector  
%          ev         element mapping matrix
%   output
%          Q          mass matrix 
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
fprintf('setting up Q1 mass matrices...  ')
%
% initialise global matrices
      r = sparse(nvtx,nvtx);
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
        re = zeros(nel,4,4);
%  loop over 2x2 Gauss points
         for igpt = 1:4
         sigpt=s(igpt);
         tigpt=t(igpt);
%  evaluate derivatives etc
         [jac,invjac,phi,dphidx,dphidy] = deriv(sigpt,tigpt,xl_v,yl_v);
            for j = 1:4
               for i = 1:4
               re(:,i,j) = re(:,i,j)  + phi(:,i).*phi(:,j) .* jac(:);
               end
	        end
% end of Gauss point loop
         end
% perform assembly of global matrix  and source vector 
      for krow=1:4
	  nrow=ev(:,krow);	 
          for kcol=1:4
		  ncol=ev(:,kcol);	  
          r = r + sparse(nrow,ncol,re(:,krow,kcol),nvtx,nvtx);
          end
      end
%
Q = r;
return

