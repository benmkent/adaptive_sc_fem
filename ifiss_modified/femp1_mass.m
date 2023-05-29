function [q] = femp1_mass(xy,evt)
%FEMP1_DIFF  set up linear anisotropic diffusion matrices
%   [A,Q,f,Ae,Qe] = femp1_diff(xy,evt);
%   input
%          xy         vertex coordinate vector  
%          evt        element mapping matrix
%   output
%          Q          mass matrix 
%
%   Natural boundary conditions apply. Dirichlet conditions
%   must be explicitly enforced by calling function nonzerobc.
%    TIFISS function: DJS; 3 March 2017.
% Copyright (c) 2007 C.E. Powell, D.J. Silvester
x=xy(:,1); y=xy(:,2);
nvtx=length(x);
nel=length(evt(:,1));
fprintf('setting up P1 mass matrices...  ')
%
% initialise global matrices
      q = sparse(nvtx,nvtx);
%
% inner loop over elements    
        for ivtx = 1:3
        xl_v(:,ivtx) = x(evt(:,ivtx));
        yl_v(:,ivtx) = y(evt(:,ivtx)); 
	    end
        qe = zeros(nel,3,3);
%
%------------------------------------------------
% 3 point Gauss rule integration
nngpt=7; [s,t,wt]=triangular_gausspoints(nngpt);

%
% loop over Gauss points
      for igpt = 1:nngpt
         sigpt=s(igpt);
         tigpt=t(igpt);
%  evaluate derivatives etc
         [jac,invjac,phi,dphidx,dphidy] = tderiv(sigpt,tigpt,xl_v,yl_v);
         for j = 1:3
               for i = 1:3
               qe(:,i,j) = qe(:,i,j)  + wt(igpt)*phi(:,i).*phi(:,j) .* jac(:);
               end
 	        end
% end of Gauss point loop
      end
%
% perform assembly of global matrix  and source vector
      for krow=1:3
	  nrow=evt(:,krow);	 
          for kcol=1:3
		  ncol=evt(:,kcol);	  
          q = q + sparse(nrow,ncol,qe(:,krow,kcol),nvtx,nvtx);
          end
      end
%
fprintf('done\n')
return

