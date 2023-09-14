function [jmp] = p1fluxjmps(p1sol,eex,xy,evt,eboundt,tve,s, diff)
%P1FLUXJMPS  corrected flux jumps for triangular P1 grid
%   [jmp] = p1fluxjmps(p1sol,eex,xy,evt,eboundt,tve,s)
%   input
%          p1sol      vertex solution vector
%          eex        element connectivity array
%          xy         vertex coordinate vector
%          evt        element mapping matrix
%          eboundt    element boundary mapping matrix
%          tve        edge location array
%          s          gaussian point in [0,1]
%   output
%          jmp        component elementwise edge flux jumps
%
%   calls functions:
%          gausspoints_oned
%          tderiv          
%          tgauss_adiff
% -------------------------------------------------------------------------
%    TIFISS function:
% Copyright (c) 2017 Alex Bespalov, Leonardo Rocchi

  x = xy(:,1);
  y = xy(:,2);
  nel = length(evt(:,1)); 

% initialise global matrices
  jmp = zeros(nel,3);         
  flux = zeros(nel,4);
  zero_v = zeros(nel,1);
  one_v = ones(nel,1);

% inner loop over elements    
  for ivtx = 1:3
      xl_v(:,ivtx) = x(evt(:,ivtx));
      yl_v(:,ivtx) = y(evt(:,ivtx)); 
      sl_v(:,ivtx) = p1sol(evt(:,ivtx)); 
  end
 
% first REFERENCE edge - diagonal
  [jac,invjac,phi,dphidx,dphidy] = tderiv(s,1-s,xl_v,yl_v);
  [diffx,diffy] = tgauss_adiff_fnhandle(s,1-s,xl_v,yl_v,diff); 
% first PHYSICAL edge
  hx_v = xl_v(:,3) - xl_v(:,2);
  hy_v = yl_v(:,3) - yl_v(:,2);
  he_v = sqrt(hx_v.*hx_v + hy_v.*hy_v);      % length of the edge
  sx_v = hx_v./he_v;     sy_v = hy_v./he_v;  % unit tangential components
  nx_v = sy_v;           ny_v = -sx_v;       % unit normal components   
  fx_v = zero_v;
  for ivtx = 1:3
      fx_v = fx_v + sl_v(:,ivtx).*( diffx(:).*dphidx(:,ivtx).*nx_v + ...
                             diffy(:).*dphidy(:,ivtx).*ny_v ).*invjac(:);
  end
  flux(:,1) = fx_v;
   

% second REFERENCE edge - left
  [jac,invjac,phi,dphidx,dphidy] = tderiv(0,s,xl_v,yl_v);
  [diffx,diffy] = tgauss_adiff_fnhandle(s,1-s,xl_v,yl_v,diff); 
% second PHYSICAL edge
  hx_v = xl_v(:,1) - xl_v(:,3);
  hy_v = yl_v(:,1) - yl_v(:,3);
  he_v = sqrt(hx_v.*hx_v + hy_v.*hy_v);
  sx_v = hx_v./he_v;     sy_v = hy_v./he_v; 
  nx_v = sy_v;           ny_v = -sx_v;
  fx_v = zero_v;
  for ivtx = 1:3
      fx_v = fx_v + sl_v(:,ivtx).*( diffx(:).*dphidx(:,ivtx).*nx_v + ...
                            diffy(:).*dphidy(:,ivtx).*ny_v ).*invjac(:);
  end
  flux(:,2) = fx_v;


% third REFERENCE edge - bottom
  [jac,invjac,phi,dphidx,dphidy] = tderiv(s,0,xl_v,yl_v);
  [diffx,diffy] = tgauss_adiff_fnhandle(s,1-s,xl_v,yl_v,diff); 
% third PHYSICAL edge
  hx_v = xl_v(:,2) - xl_v(:,1); 
  hy_v = yl_v(:,2) - yl_v(:,1);
  he_v = sqrt(hx_v.*hx_v + hy_v.*hy_v); 
  sx_v = hx_v./he_v;     sy_v = hy_v./he_v;  
  nx_v = sy_v;           ny_v = -sx_v; 
  fx_v = zero_v;
  for ivtx = 1:3
      fx_v = fx_v + sl_v(:,ivtx).*( diffx(:).*dphidx(:,ivtx).*nx_v + ...
                             diffy(:).*dphidy(:,ivtx).*ny_v ).*invjac(:);
  end
  flux(:,3) = fx_v;
       
% add zero column for boundary jumps
  flux(:,4) = zero_v;

% replace zero indices in array tve by 4s
  tvx = tve;
  tvx(find(tve==0))=4;
      
% evaluate flux jump on each edge in turn
% A(sub2ind(size(A),ii,jj)) pulls out the entries of flux indexed by ii and jj

% first edge
  jmp(:,1) = flux(:,1) + flux( sub2ind([nel,4],eex(:,1),tvx(:,1)) );
% second edge
  jmp(:,2) = flux(:,2) + flux( sub2ind([nel,4],eex(:,2),tvx(:,2)) );
% third edge
  jmp(:,3) = flux(:,3) + flux( sub2ind([nel,4],eex(:,3),tvx(:,3)) );
   
% remove Dirichlet boundary edge contributions
  nbde = length(eboundt(:,1));
  for k=1:nbde
      jmp(eboundt(k,1),eboundt(k,2))=0;
  end

end  % end function
