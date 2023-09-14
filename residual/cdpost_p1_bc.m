function [err_p,elerr_p] = cdpost_p1_bc(aez,fez,elerror,xy,evt,eboundt,bc_fn)
%DIFFPOST_P1_BC postprocesses Poisson error estimator at boundary elements
% [err_p,elerr_p] = diffpost_p1_bc(aez,fez,elerror,xy,evt,eboundt);
%   input
%          aez       elementwise Poisson problem matrices
%          fez       elementwise rhs vectors
%          elerror   elementwise error estimate (without BC imposition) 
%          xy        vertex coordinate vector  
%          ev        element mapping matrix
%          eboundt   element edge boundary matrix 
%   output
%          err_p     global error estimate 
%          elerr_p   elementwise error estimate
% -------------------------------------------------------------------------
%    TIFISS function:
% Copyright (c) 2017 Alex Bespalov, Leonardo Rocchi

% NOTE: modification
% -------------------------------------------------------------------------
% The only difference with the original diffpost_p1_bc.m is in the command
% 
%    ae = squeeze(aez(el1e,1:4,1:4));
% 
% that here becomes
%
%    ae = squeeze(aez(el1e,:,:));
%
% without any specification on numbers.

  x = xy(:,1); 
  y = xy(:,2);
  nel = size(evt,1);
  lev = [evt,evt(:,1),evt(:,2)];  
  elerr_p = elerror.*elerror;
  
% Recompute contributions from elements with Dirichlet boundaries
  nbde = size(eboundt,1);
  ebdy = zeros(nel,1);
  edge = zeros(nel,1);
% isolate boundary elements
  for el = 1:nbde
      ee = eboundt(el,1);
      ebdy(ee) = ebdy(ee)+1; 
      edge(ee) = eboundt(el,2);
  end  
         
% Two edge elements
  k2 = find(ebdy==2);
  nel2b = length(k2);
  % loop over two edge elements
  for el = 1:nel2b
      el2e = k2(el);
      kk = find(eboundt(:,1) == el2e);
      edges = eboundt(kk,2);
      % set up original matrix and RHS vector
      ae = squeeze(aez(el2e,:,:)); 
      fe = fez(el2e,:)';
      % set up local coordinates and impose interpolated error as Dirichlet bc
      xl = x(lev(el2e,:)); 
      yl = y(lev(el2e,:)); 
      [bae,fe] = localbc_p_bc(ae,fe,edges,xl,yl,bc_fn);
      fprintf('\n<strong>Warning:</strong> element %g has two boundary edges\n',el2e)
      % solve local problem
      err = bae\fe;
      elerr_p(el2e,1) = err'*fe;
  end
% end of element loop

% One edge elements
  k1 = find(ebdy==1);
  nel1b = length(k1);
% loop over one edge elements
  for el = 1:nel1b       
      el1e = k1(el);
      kk = find(eboundt(:,1) == el1e);
      edges = eboundt(kk,2);
      % set up original matrix and RHS vector 
      fe = fez(el1e,:)';
	  ae = squeeze(aez(el1e,:,:));
      % set up local coordinates and impose interpolated error as Dirichlet bc
      xl = x(lev(el1e,:));
      yl = y(lev(el1e,:));
      [bae,fe] = localbc_p_bc(ae,fe,edges,xl,yl,bc_fn);
      % solve local problem
      err = bae\fe;
      elerr_p(el1e,1) = err'*fe;
  end
% end of element loop
     
% Final error
  err_p = sqrt(sum(elerr_p));
  elerr_p = sqrt(elerr_p);
  fprintf('boundary correction done\n');
         
end  % end function
