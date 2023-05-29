function [a,n,r,epe,eph,epw, n2, xl_v, yl_v,flowx, flowy] = femq1_cd_bk(xy,ev, wind_fn,varargin)
%FEMQ1_CD vectorized bilinear coefficient matrix generator
%    [A,N,Q,epe,eph,epw] = femq1_cd(xy,ev);
%   input
%          xy       vertex coordinate vector  
%          ev       element mapping matrix
%   output
%          A        stiffness matrix
%          N        convection matrix
%          Q        mass matrix 
%          epe      viscosity normalised element peclet numbers 
%          eph      flow specific element lengths 
%          epw      centroid evaluated wind 
%
%   Natural boundary conditions apply. Dirichlet conditions
%   must be explicitly enforced by calling function nonzerobc.
%   IFISS function: DJS; 5 March 2005.
% Copyright (c) 2005 D.J. Silvester, H.C. Elman, A. Ramage

if nargin == 3
    diff_fn = @(x1,x2) ones(size(x1,1),1);
else
    diff_fn = varargin{1};
end

x=xy(:,1); y=xy(:,2);
nvtx=length(x);
nel=length(ev(:,1));
lx=max(x)-min(x); ly=max(y)-min(y);
hx=max(diff(x)); hy=max(diff(y));
fprintf('setting up Q1 convection-diffusion matrices...  ')
%
% initialise global matrices
      a = sparse(nvtx,nvtx);
      n = sparse(nvtx,nvtx);
      n2 = sparse(nvtx,nvtx);
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
      ae = zeros(nel,4,4);
	  ne = zeros(nel,4,4);
	  n2e = zeros(nel,4,4);
      re = zeros(nel,4,4);
      fe = zeros(nel,4);  
% loop over 2x2 Gauss points
      for igpt = 1:4
         sigpt=s(igpt);
         tigpt=t(igpt);
%  evaluate derivatives etc
         [jac,invjac,phi,dphidx,dphidy] = deriv(sigpt,tigpt,xl_v,yl_v);
         [flowx,flowy] = gauss_transprt_modified(sigpt,tigpt,xl_v,yl_v,wind_fn);
         diffeval = diff_eval_modified(sigpt, tigpt, xl_v, yl_v, diff_fn);
         for j = 1:4
            for i = 1:4
               ae(:,i,j) = ae(:,i,j)  + diffeval.*dphidx(:,i).*dphidx(:,j) .* invjac(:);
               ae(:,i,j) = ae(:,i,j)  + diffeval.*dphidy(:,i).*dphidy(:,j) .* invjac(:);
               re(:,i,j) = re(:,i,j)  + phi(:,i).*phi(:,j) .* jac(:);
               ne(:,i,j) = ne(:,i,j) + flowx(:) .* phi(:,i) .* dphidx(:,j);
               ne(:,i,j) = ne(:,i,j) + flowy(:) .* phi(:,i) .* dphidy(:,j);
               n2e(:,i,j) = n2e(:,i,j) + flowx(:) .* dphidx(:,i) .* flowx(:) .* dphidx(:,j);
               n2e(:,i,j) = n2e(:,i,j) + flowy(:) .* dphidy(:,i) .* flowy(:) .* dphidy(:,j);
               n2e(:,i,j) = n2e(:,i,j) + 2 .* flowx(:) .* dphidx(:,i) .* flowy(:) .* dphidy(:,j);
            end
	     end
% end of Gauss point loop
      end
%
% perform assembly of global matrix  and source vector 
      for krow=1:4
	     nrow=ev(:,krow);	 
         for kcol=1:4
            ncol=ev(:,kcol);	  
            a = a + sparse(nrow,ncol,ae(:,krow,kcol),nvtx,nvtx);
            r = r + sparse(nrow,ncol,re(:,krow,kcol),nvtx,nvtx);
            n = n + sparse(nrow,ncol,ne(:,krow,kcol),nvtx,nvtx);
            n2 = n2 + sparse(nrow,ncol,n2e(:,krow,kcol),nvtx,nvtx);
         end
      end
%
% computation of element Peclet number (at the centroid)         
% rectangle specific calculation here
      hx=abs(xl_v(:,2)-xl_v(:,1)); hy=abs(yl_v(:,3)-yl_v(:,2));
      [flowx,flowy] = gauss_transprt_modified(0,0,xl_v,yl_v,wind_fn);
      flow_l2 = sqrt(flowx(:) .* flowx(:) + flowy(:) .* flowy(:));
      if     all(flowx==0), flow_h=hy;
	  elseif all(flowy==0), flow_h=hx;
      else
         angle = atan(abs(flowy./flowx));
         flow_h = min([hx./cos(angle),hy./sin(angle)],[],2);
      end
      eph = flow_h;
      epe = flow_h.*flow_l2/2;
	  epw = flow_l2;
%
fprintf('done\n')
return
