function [elerr_p,fe,AE] = diffpost_p1_with_p2_for_wind(xy,evt,eex,tve,els,eboundt,p1sol, dU, w_fn)
%DIFFPOST_P1_WITH_P2 a posteriori estimation for P1 using P2 bubble functions
%
%   [elerr_p,fe,AE] = diffpost_p1_with_p2(xy,evt,eex,tve,els,eboundt,p1sol)
%
%   input:
%               xy    vertex coordinate vector  
%              evt    element mapping matrix
%              eex    element connectivity array
%              tve    edge location array
%              els    elementwise edge lengths
%          eboundt    element boundary mapping matrix
%            p1sol    P1 solution for diffusion problem
%
%   output:
%          elerr_p    elementwise error estimate
%               fe    elementwise rhs vectors
%               AE    elementwise Poisson problem matrices
%
% Function(s) called: triangular_gausspoints
%                     tderiv
%                     tqderiv
%                     tgauss_adiff
%                     intres_p1_with_p2
%                     edgeres_p1_with_p2
%
% See also DIFFPOST_P1_WITH_P1
%
% Last update: 01/02/2017
% -------------------------------------------------------------------------
%    TIFISS function:
% Copyright (c) 2017 Alex Bespalov, Leonardo Rocchi

  fprintf('P1 local error estimator using bubble functions...\n');
  x = xy(:,1);
  y = xy(:,2);
  nel = length(evt(:,1));
  elerr_p = zeros(nel,1);

% Number of bubble functions  
  nnode = 4; 
             
% Construct 2D gaussian rule over the reference triangle
  nngpt = 7;
  [s,t,wt] = triangular_gausspoints(nngpt);

  xl_v = zeros(nel,3);
  yl_v = zeros(nel,3);
% Recover local coordinates  
  for ivtx = 1:3 
      xl_v(:,ivtx) = x(evt(:,ivtx));
      yl_v(:,ivtx) = y(evt(:,ivtx));
  end

% Initialisation 
  ae = zeros(nel,nnode,nnode);

% LHS of the linear system
% -------------------------------------------------------------------------

% Loop over Gauss points
  for igpt = 1:nngpt
      sigpt = s(igpt); 
      tigpt = t(igpt);
      wght = wt(igpt);  
      % Evaluate derivatives
      [~,invjac_v,~,~,~] = tderiv(sigpt,tigpt,xl_v,yl_v);
      [~,dpsidx_v,dpsidy_v] = tqderiv(sigpt,tigpt,xl_v,yl_v);
%       % Evaluate variable diffusion coefficients
%       [diffx,diffy] = tgauss_adiff(sigpt,tigpt,xl_v,yl_v);
      % Loop over the four bubble functions
      for j = 1:nnode
          for i = 1:nnode
              ae(:,i,j) = ae(:,i,j) + wght .* dpsidx_v(:,i+3) .* dpsidx_v(:,j+3) .* invjac_v(:);
              ae(:,i,j) = ae(:,i,j) + wght .* dpsidy_v(:,i+3) .* dpsidy_v(:,j+3) .* invjac_v(:);
          end
      end
      % end bubble functions loop
  end
% end of Gauss point loop


% Saving the LHS matrix of the system for output before factorization
  AE = ae;
   
% RHS of the linear system
% -------------------------------------------------------------------------

% Element residual
% This assumes diffusion coefficent has zero gradient currently
  [res_int] = intres_p1_with_p2_windforcing(xy,evt,p1sol, dU, w_fn);
    
% Edge residual
  [res_edge] = edgeres_p1_with_p2(xy,evt,eboundt,p1sol,eex,tve,els);
  
  fprintf('internal_res = %7.4e;     edge_res = %7.4e\n',norm(res_int), norm(res_edge));
  
% Final rhs of the linear system
  for j = 1:3
      res_int(:,j) = res_int(:,j) - res_edge(:,j);
  end
  fe = res_int; 

% Vectorized code - LDLT factorization
  nn = nnode;
  dd = zeros(nel,nn); 
  rr = zeros(nel,nn);

  for kk=1:nn-1  
      for pp = 1:kk-1;
          rr(1:nel,pp) = dd(1:nel,pp).*ae(1:nel,kk,pp);
      end
      dd(1:nel,kk) = ae(1:nel,kk,kk);
      for pp = 1:kk-1;
          dd(1:nel,kk)= dd(1:nel,kk) - ae(1:nel,kk,pp).*rr(1:nel,pp);
      end
      for ii = kk+1:nn
          for pp = 1:kk-1;
              ae(1:nel,ii,kk) = ae(1:nel,ii,kk) - ae(1:nel,ii,pp).*rr(1:nel,pp);
          end
         ae(1:nel,ii,kk) = ae(1:nel,ii,kk)./dd(1:nel,kk);
      end 
  end
  
  for pp = 1:nn-1
      rr(1:nel,pp) = dd(1:nel,pp).*ae(1:nel,nn,pp);
  end

  dd(1:nel,nn) = ae(1:nel,nn,nn);

  for pp = 1:nn-1;
      dd(1:nel,nn) = dd(1:nel,nn)- ae(1:nel,nn,pp).*rr(1:nel,pp);
  end

% overwrite diagonal entries
  for kk=1:nn
      ae(1:nel,kk,kk) = dd(1:nel,kk);
  end

% forward-backward substitutions ...
  xx = element_lusolve(ae,fe);
  elerr = xx';
                   
  for ivtx = 1:nnode
      elerr_p(:) = elerr_p(:) + fe(:,ivtx) .* elerr(ivtx,:)';
  end
  elerr_p = sqrt(elerr_p);

  fprintf('estimated energy error is %10.4e \n',norm(elerr_p,2));

end  % end function
