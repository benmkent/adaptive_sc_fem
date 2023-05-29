function n = femq2_convection(xy, ev, xflow, yflow)
%FEMQ2_CONVECTION vectorized biquadratic coefficient matrix generator
%   N = femq2_convection(xy,mv,flow,tout);
%   input
%          xy       vertex coordinate vector  
%          mv       element mapping matrix
%        flow       streamfunction defined on the Q1 grid
%        tout       output on/off switch (optional)
%   output
%          N        Q2 convection matrix
%
%   Natural boundary conditions apply. Dirichlet conditions
%   must be explicitly enforced by calling function nonzerobc.
%   s-IFISS function: DJS; 3 November 2014.
% Copyright (c) 2014 D.J. Silvester, H.C. Elman
if nargin < 4, tout = 1; end
x=xy(:,1); y=xy(:,2); 
nvtx=length(x); nel=length(ev(:,1));
if tout==1,
fprintf('assembling Q1 convection matrix ...\n')
end

% initialise global matrices
      n = sparse(nvtx,nvtx);
%
% % set up 2x2 Gauss points
%     nngpt=4;
% 
%       wt(1:4)=1;
%       gpt=1.0e0/sqrt(3.0e0);
%       s(1) = -gpt;  t(1) = -gpt;
%       s(2) =  gpt;  t(2) = -gpt;
%       s(3) =  gpt;  t(3) =  gpt;
%       s(4) = -gpt;  t(4) =  gpt;

% set up 3x3 Gauss points
nngpt=9;

gpt=sqrt(0.6);
s(1) = -gpt; t(1) = -gpt; wt(1)=25/81;
s(2) =  gpt; t(2) = -gpt; wt(2)=25/81;
s(3) =  gpt; t(3) =  gpt; wt(3)=25/81;
s(4) = -gpt; t(4) =  gpt; wt(4)=25/81;
s(5) =  0.0; t(5) = -gpt; wt(5)=40/81;
s(6) =  gpt; t(6) =  0.0; wt(6)=40/81;
s(7) =  0.0; t(7) =  gpt; wt(7)=40/81;
s(8) = -gpt; t(8) =  0.0; wt(8)=40/81;
s(9) =  0.0; t(9) =  0.0; wt(9)=64/81;

% inner loop over elements    
      for ivtx = 1:4
         xl_v(:,ivtx) = x(ev(:,ivtx));
         yl_v(:,ivtx) = y(ev(:,ivtx));
         fl_v(:,ivtx) = flow(ev(:,ivtx));
      end
	  ne = zeros(nel,4,4);
% loop over 2x2 Gauss points
      for igpt = 1:9
         sigpt=s(igpt);
         tigpt=t(igpt);
%  evaluate derivatives etc
         [jac,invjac,phi,dphidx,dphidy] = deriv(sigpt,tigpt,xl_v,yl_v);
         [psi,dpsidx,dpsidy] = qderiv(sigpt,tigpt,xl_v,yl_v); 

%           flowx = zeros(nel,1); flowy=zeros(nel,1);
%           for k=1:4 
% %           flowx(:) = flowx(:) - fl_v(:,k) .* dphidy(:,k);
% %           flowy(:) = flowy(:) + fl_v(:,k) .* dphidx(:,k);
% 
%           flowx(:) = flowx(:) - fl_v(:,k) .* dpsidx(:,k);
%           flowy(:) = flowy(:) + fl_v(:,k) .* dpsidy(:,k);
%           end

          flowx = xflow()

          for j = 1:4
              for i = 1:4
               ne(:,i,j) = ne(:,i,j) + flowx(:) .* dphidx(:,j)  .* phi(:,i);
               ne(:,i,j) = ne(:,i,j) + flowy(:)  .* dphidy(:,j) .* phi(:,i);
              end
          end
% end of Gauss point loop
      end
%
%%  element assembly into global matrix
      for krow=1:4
      nrow=ev(:,krow);
          for kcol=1:4
          ncol=ev(:,kcol);
            n = n + sparse(nrow,ncol,ne(:,krow,kcol),nvtx,nvtx);
          end
      end
      n=(n-n')/2;
%
%
%if tout==1, fprintf('done.\n'), end
return
